#include <doctest/doctest.h>

#include "Wave.hpp"

#include <type_traits>
#include <utility>

// Wave's copy operations are shallow: they share the audio buffer and copy the metadata
// and content revision (see the class docs in Wave.hpp). Revision behaviour is covered
// in WaveRevisionTests.cpp; these cases pin the operator signatures and buffer sharing.

static_assert(std::is_same_v<decltype(std::declval<wm::Wave&>() = std::declval<const wm::Wave&>()),
                             wm::Wave&>,
              "copy assignment must return Wave& so assignments can be chained");

TEST_CASE("copy assignment returns a reference to the assigned-to wave")
{
    wm::Wave source = wm::Wave::generateRandom(100);
    wm::Wave target;
    wm::Wave& result = (target = source);
    CHECK(&result == &target);
}

TEST_CASE("chained copy assignment shares one buffer and format across all waves")
{
    wm::Wave source(64, 1, 24, 48000);
    wm::Wave middle = wm::Wave::generateRandom(100);
    wm::Wave last;

    last = middle = source;

    CHECK(middle.constAudioData() == source.constAudioData());
    CHECK(last.constAudioData() == source.constAudioData());
    CHECK(last.getNumFrames() == 64);
    CHECK(last.getNumChannels() == 1);
    CHECK(last.getSampleBitDepth() == 24);
    CHECK(last.getSampleRate() == 48000);
    CHECK(last.contentRevision() == source.contentRevision());
}

TEST_CASE("self-assignment leaves the wave unchanged")
{
    wm::Wave wave = wm::Wave::generateRandom(100);
    const float* data = wave.constAudioData();
    const uint32_t frames = wave.getNumFrames();
    const uint64_t revision = wave.contentRevision();

    wm::Wave& alias = wave;
    wave = alias;

    CHECK(wave.constAudioData() == data);
    CHECK(wave.getNumFrames() == frames);
    CHECK(wave.contentRevision() == revision);
}

TEST_CASE("copy construction shares the source buffer")
{
    wm::Wave source = wm::Wave::generateRandom(100);
    wm::Wave copy(source);
    CHECK(copy.constAudioData() == source.constAudioData());
    CHECK(copy.getNumFrames() == source.getNumFrames());
}
