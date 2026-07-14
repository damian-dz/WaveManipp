#include <doctest/doctest.h>

#include "Dsp.hpp"
#include "Wave.hpp"

#include <algorithm>
#include <vector>

namespace {

wm::Wave makeTestWave(uint32_t frames, uint16_t channels, uint32_t sampleRate)
{
    wm::Wave wave(frames, channels, 16, sampleRate);
    std::vector<float> samples(static_cast<size_t>(frames) * channels);
    for (uint32_t frame = 0; frame < frames; ++frame) {
        for (uint16_t channel = 0; channel < channels; ++channel) {
            const float impulse =
                (frame == 120 + channel * 17) ? 0.8f : 0.f;
            const float repeating =
                static_cast<float>((frame * 13 + channel * 7) % 29 - 14) * 0.01f;
            samples[static_cast<size_t>(frame) * channels + channel] =
                impulse + repeating;
        }
    }
    wave.setAudio(samples.data(), static_cast<uint32_t>(samples.size()));
    return wave;
}

std::vector<float> copySamples(const wm::Wave& wave)
{
    const float* data = wave.constAudioData();
    return std::vector<float>(data, data + wave.getNumSamples());
}

} // namespace

TEST_CASE("ReverbProcessor matches range reverb across arbitrary block boundaries")
{
    constexpr uint32_t frames = 2600;
    constexpr uint16_t channels = 2;
    constexpr uint32_t sampleRate = 48000;
    constexpr uint32_t startFrame = 95;
    constexpr uint32_t endFrame = 2300;

    wm::dsp::ReverbParams params;
    params.roomSize = 0.73f;
    params.damping = 0.35f;
    params.wetMix = 0.42f;

    wm::Wave expected = makeTestWave(frames, channels, sampleRate);
    wm::Wave source(expected);
    wm::dsp::reverb(expected, startFrame, endFrame, params);

    std::vector<float> actual = copySamples(source);
    wm::dsp::ReverbProcessor processor(sampleRate, channels, params);
    uint32_t offset = startFrame;
    const uint32_t blockSizes[] = { 1, 37, 128, 13, 251 };
    size_t blockIndex = 0;
    while (offset < endFrame) {
        const uint32_t requested = blockSizes[blockIndex++ % 5];
        const uint32_t blockFrames = std::min(requested, endFrame - offset);
        processor.process(actual.data() + static_cast<size_t>(offset) * channels,
                          blockFrames);
        offset += blockFrames;
    }

    const float* expectedSamples = expected.constAudioData();
    REQUIRE(actual.size() == expected.getNumSamples());
    for (size_t sample = 0; sample < actual.size(); ++sample) {
        CHECK(actual[sample] == doctest::Approx(expectedSamples[sample]).epsilon(0.000001));
    }
}
