#include <doctest/doctest.h>

#include "Dsp.hpp"
#include "Wave.hpp"

#include <filesystem>
#include <vector>

// Phase 2 of docs/WAVEFORM_PEAK_CACHE_PHASE2_PLAN.md section 1: wm::Wave::contentRevision()
// is the explicit source-invalidation contract the background peak-build queue relies on.
// These tests exercise WaveManipp directly (no TrackCtrl/AudioClip involvement), which is
// why they run from the WaveManipp_tests target rather than the app's WaveformPeakCache_tests.

namespace {

// Brackets an operation with two default-constructed marker waves to count exactly how
// many global revision ticks the operation itself consumed, independent of whatever
// absolute revision value the shared atomic counter happens to be at (other tests share
// it too). Markers must be default-constructed specifically: the numFrames constructor
// (as used by e.g. generateRandom()) burns three ticks per call (its own init-list, plus
// reserveMemory(), plus setSampleRate() -- all individually correct, just not minimal),
// which would silently break the "credited ticks" arithmetic below. before.contentRevision()
// and after.contentRevision() straddle the operation; the operation itself is credited with
// (after - before - 1) ticks, since `after`'s own construction consumes exactly one more
// tick than whatever came right before it.
uint64_t ticksConsumedBetween(const wm::Wave& before, const wm::Wave& after)
{
    return after.contentRevision() - before.contentRevision() - 1;
}

} // namespace

// ============================================================================
// Fresh content and shallow-copy preservation
// ============================================================================

TEST_CASE("new waves receive distinct revisions")
{
    wm::Wave a = wm::Wave::generateRandom(100);
    wm::Wave b = wm::Wave::generateRandom(100);
    CHECK(a.contentRevision() != b.contentRevision());
}

TEST_CASE("independently constructed waves with identical PCM are not required to share a revision")
{
    wm::Wave a(10, 1);
    wm::Wave b(10, 1); // both default-initialized (silent) -- identical content
    CHECK(a.contentRevision() != b.contentRevision());
}

TEST_CASE("copy construction preserves the source's revision")
{
    wm::Wave original = wm::Wave::generateRandom(100);
    const uint64_t revision = original.contentRevision();
    wm::Wave copy(original);
    CHECK(copy.contentRevision() == revision);
}

TEST_CASE("operator= preserves the assigned-from revision")
{
    wm::Wave source = wm::Wave::generateRandom(100);
    wm::Wave target = wm::Wave::generateRandom(50);
    REQUIRE(target.contentRevision() != source.contentRevision());
    target = source;
    CHECK(target.contentRevision() == source.contentRevision());
}

// ============================================================================
// Every mutating call bumps the revision exactly once
// ============================================================================

TEST_CASE("audioData() changes the revision")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    wave.audioData()[0] = 0.5f;
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("changeVolume() changes the revision")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    wave.changeVolume(0.5f);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("reverse() changes the revision")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    wave.reverse();
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("insertAudio() changes the revision")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    std::vector<float> patch(wave.getNumSamples(), 0.25f);
    wave.insertAudio(0, patch);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("setAudio() changes the revision -- same size (in-place) branch")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    std::vector<float> replacement(wave.getNumSamples(), 0.1f);
    wave.setAudio(replacement);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("setAudio() changes the revision -- different size (resize) branch")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    std::vector<float> replacement(50, 0.1f);
    wave.setAudio(replacement);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("swapChannels() changes the revision")
{
    wm::Wave wave(100, 2);
    const uint64_t before = wave.contentRevision();
    wave.swapChannels(0, 1);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("downmixToMono() changes the revision")
{
    wm::Wave wave(100, 2);
    const uint64_t before = wave.contentRevision();
    wave.downmixToMono();
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("upmixToStereo() changes the revision")
{
    wm::Wave wave(100, 1);
    const uint64_t before = wave.contentRevision();
    wave.upmixToStereo();
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("append() (resizeMemory) changes the revision")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    wm::Wave other = wm::Wave::generateRandom(50);
    const uint64_t before = wave.contentRevision();
    wave.append(other);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("setSampleRate() changes the revision")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t before = wave.contentRevision();
    wave.setSampleRate(wave.getSampleRate() + 1);
    CHECK(wave.contentRevision() != before);
}

TEST_CASE("setSampleBitDepth() does not change the revision")
{
    wm::Wave wave(100, 1, 16);
    const uint64_t before = wave.contentRevision();
    wave.setSampleBitDepth(24);
    CHECK(wave.contentRevision() == before);
}

TEST_CASE("loading the same file twice produces two distinct revisions")
{
    wm::Wave source = wm::Wave::generateRandom(200);
    const auto path = std::filesystem::temp_directory_path() / "wavemanipp_revision_test.wav";
    source.saveAs(path.string());

    wm::Wave loadedA(path.string());
    wm::Wave loadedB(path.string());
    CHECK(loadedA.contentRevision() != loadedB.contentRevision());

    std::filesystem::remove(path);
}

// ============================================================================
// Per-sample DSP operations bump the revision exactly once, not once per sample
// ============================================================================

TEST_CASE("applyGainRange bumps the revision exactly once across many samples")
{
    wm::Wave wave = wm::Wave::generateRandom(50000); // large enough that a per-sample bug is obvious
    wm::Wave before; // default-constructed: the only marker guaranteed to consume exactly one tick
    wm::dsp::applyGainRange(wave, 0, wave.getNumFrames(), 1.5f);
    wm::Wave after; // likewise
    CHECK(ticksConsumedBetween(before, after) == 1);
}

TEST_CASE("fadeIn bumps the revision exactly once across many samples")
{
    wm::Wave wave = wm::Wave::generateRandom(50000);
    wm::Wave before; // default-constructed: the only marker guaranteed to consume exactly one tick
    wm::dsp::fadeIn(wave, 0, wave.getNumFrames(), wm::dsp::FadeCurve::Linear);
    wm::Wave after; // likewise
    CHECK(ticksConsumedBetween(before, after) == 1);
}

TEST_CASE("processBiquadRange bumps the revision exactly once across many samples")
{
    wm::Wave wave = wm::Wave::generateRandom(50000);
    const auto coeffs = wm::dsp::makePeaking(1000.f, static_cast<float>(wave.getSampleRate()), 3.f, 0.7f);
    wm::Wave before; // default-constructed: the only marker guaranteed to consume exactly one tick
    wm::dsp::processBiquadRange(wave, 0, wave.getNumFrames(), coeffs);
    wm::Wave after; // likewise
    CHECK(ticksConsumedBetween(before, after) == 1);
}

TEST_CASE("reverseChannel bumps the revision exactly once across many samples")
{
    wm::Wave wave = wm::Wave::generateRandom(50000);
    wm::Wave before; // default-constructed: the only marker guaranteed to consume exactly one tick
    wm::dsp::reverseChannel(wave, 0);
    wm::Wave after; // likewise
    CHECK(ticksConsumedBetween(before, after) == 1);
}

// ============================================================================
// COW sharing and undo-style revision recovery
// ============================================================================

TEST_CASE("mutating one COW copy changes only that copy's revision")
{
    wm::Wave original = wm::Wave::generateRandom(100);
    wm::Wave sibling(original); // shares the buffer, shares the revision
    const uint64_t sharedRevision = original.contentRevision();
    REQUIRE(sibling.contentRevision() == sharedRevision);

    original.changeVolume(0.5f); // detaches `original` onto a new buffer

    CHECK(original.contentRevision() != sharedRevision);
    CHECK(sibling.contentRevision() == sharedRevision); // untouched
}

TEST_CASE("edit, undo via reassignment, then a different edit produces three distinct revisions")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const uint64_t originalRevision = wave.contentRevision();
    wm::Wave snapshot(wave); // undo snapshot, shares the original revision

    wave.changeVolume(0.5f); // "edit"
    const uint64_t editedRevision = wave.contentRevision();
    CHECK(editedRevision != originalRevision);

    wave = snapshot; // "undo" -- restores the original revision, not a fresh one
    CHECK(wave.contentRevision() == originalRevision);

    wave.reverse(); // "a different edit"
    const uint64_t secondEditRevision = wave.contentRevision();
    CHECK(secondEditRevision != originalRevision);
    CHECK(secondEditRevision != editedRevision);
}
