#include <doctest/doctest.h>

#include <Dsp.hpp>

#include <vector>

namespace {

wm::Wave makeWave(uint32_t frames,
                  uint16_t channels,
                  uint32_t sampleRate,
                  const std::vector<float>& samples)
{
    wm::Wave wave(frames, channels, 32, sampleRate);
    wave.setAudio(samples.data(), static_cast<uint32_t>(samples.size()));
    return wave;
}

void checkSamplesApproxEqual(const std::vector<float>& actual,
                             const wm::Wave& expected)
{
    REQUIRE(actual.size() ==
            static_cast<std::size_t>(expected.getNumFrames()) *
                expected.getNumChannels());
    const float* expectedSamples = expected.constAudioData();
    for (std::size_t sample = 0; sample < actual.size(); ++sample) {
        CHECK(actual[sample] ==
              doctest::Approx(expectedSamples[sample]).epsilon(0.000001));
    }
}

} // namespace

TEST_CASE("PhaseVocoderTimeStretcher matches stretchTempo for mono compression")
{
    const std::vector<float> samples{
        0.00f, 0.25f, 0.50f, 0.75f,
        1.00f, 0.60f, 0.20f, -0.20f,
        -0.60f, -1.00f, -0.50f, 0.00f
    };
    wm::Wave expected = makeWave(12, 1, 44100, samples);
    wm::dsp::stretchTempo(expected, 0.75f, 8);

    wm::dsp::PhaseVocoderTimeStretcher stretcher(8);
    const std::vector<float> actual =
        stretcher.render(samples.data(), 12, 1, 0.75f);

    CHECK(expected.getNumFrames() == 9);
    checkSamplesApproxEqual(actual, expected);
}

TEST_CASE("PhaseVocoderTimeStretcher matches stretchTempo for stereo expansion")
{
    const std::vector<float> samples{
        0.0f,  1.0f,
        0.2f,  0.8f,
        0.4f,  0.6f,
        0.6f,  0.4f,
        0.8f,  0.2f,
        1.0f,  0.0f,
        0.5f, -0.5f,
        0.0f, -1.0f
    };
    wm::Wave expected = makeWave(8, 2, 48000, samples);
    wm::dsp::stretchTempo(expected, 1.5f, 8);

    wm::dsp::PhaseVocoderTimeStretcher stretcher(8);
    const std::vector<float> actual =
        stretcher.render(samples.data(), 8, 2, 1.5f);

    CHECK(expected.getNumFrames() == 12);
    checkSamplesApproxEqual(actual, expected);
}

TEST_CASE("PhaseVocoderTimeStretcher exposes output frame planning")
{
    wm::dsp::PhaseVocoderTimeStretcher stretcher(16);

    CHECK(stretcher.fftSize() == 16);
    CHECK(stretcher.outputFrameCount(10, 0.5f) == 5);
    CHECK(stretcher.outputFrameCount(10, 1.5f) == 15);
    CHECK(stretcher.outputFrameCount(1, 0.1f) == 1);
    CHECK(stretcher.outputFrameCount(0, 2.f) == 0);
}

TEST_CASE("PhaseVocoderTimeStretcher returns input for no-op factors")
{
    const std::vector<float> samples{0.f, 1.f, 2.f, 3.f};
    wm::dsp::PhaseVocoderTimeStretcher stretcher(8);

    CHECK(stretcher.render(samples.data(), 4, 1, 1.f) == samples);
    CHECK(stretcher.render(samples.data(), 4, 1, 0.f) == samples);
}
