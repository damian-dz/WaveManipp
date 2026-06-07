#include "Dsp.hpp"
#include "Fft.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

namespace wm::dsp {

namespace {

static constexpr double k_2pi = 6.283185307179586476925286766559;

// Wraps x into [-pi, pi] using the IEEE remainder.
static inline double wrapPhase(double x)
{
    return std::remainder(x, k_2pi);
}

} // anonymous namespace

void stretchTempo(Wave& wave, float factor, size_t fftSize)
{
    if (factor <= 0.f || wave.isEmpty()) return;
    if (std::abs(factor - 1.f) < 1e-5f) return;

    const size_t   N           = fftSize;
    const size_t   Ha          = N / 4;   // analysis hop  (75 % overlap)
    const size_t   Hs          = std::max<size_t>(1,
        static_cast<size_t>(std::round(Ha * static_cast<double>(factor))));

    const uint32_t inputFrames  = wave.getNumFrames();
    const int      numChannels  = static_cast<int>(wave.getNumChannels());
    const uint32_t outputFrames = static_cast<uint32_t>(std::max<int64_t>(1,
        static_cast<int64_t>(std::round(inputFrames * static_cast<double>(factor)))));

    // Number of analysis frames needed to cover all input samples.
    const size_t numAFrames = (static_cast<size_t>(inputFrames) + Ha - 1) / Ha + 2;

    // Output buffer: must hold all synthesis frames plus one extra window.
    const size_t outBufLen = numAFrames * Hs + N;

    std::vector<float>  outBuf(outBufLen * numChannels, 0.f);
    std::vector<double> normBuf(outBufLen, 0.0);

    const std::vector<float> win = makeWindow(N, WindowFunction::Hann);

    // Pre-compute normalisation: sum of squared synthesis windows at each output position.
    // This handles edge effects at the start and end automatically.
    {
        size_t pos_s = 0;
        for (size_t f = 0; f < numAFrames; ++f, pos_s += Hs) {
            for (size_t i = 0; i < N && pos_s + i < outBufLen; ++i)
                normBuf[pos_s + i] += static_cast<double>(win[i]) * win[i];
        }
    }

    const float*  src  = wave.constAudioData();
    fft::Plan     plan(N);

    std::vector<double> re(N), im(N);

    for (int ch = 0; ch < numChannels; ++ch) {
        std::vector<double> prevPhase(N / 2 + 1, 0.0);
        std::vector<double> phaseAccum(N / 2 + 1, 0.0);

        size_t pos_s = 0;
        for (size_t f = 0; f < numAFrames; ++f, pos_s += Hs) {
            const auto pos_a = static_cast<ptrdiff_t>(f * Ha);

            // --- Analysis: extract windowed frame, zero-padding at boundaries ---
            for (size_t i = 0; i < N; ++i) {
                const ptrdiff_t s = pos_a + static_cast<ptrdiff_t>(i);
                re[i] = (s >= 0 && s < static_cast<ptrdiff_t>(inputFrames))
                        ? static_cast<double>(src[s * numChannels + ch]) * win[i]
                        : 0.0;
                im[i] = 0.0;
            }

            plan.transform(re, im);

            // --- Phase vocoder: accumulate instantaneous frequencies ---
            const size_t numBins = N / 2 + 1;

            if (f == 0) {
                // Initialise synthesis phase from the first analysis frame
                // so the first output frame is phase-coherent with the input.
                for (size_t k = 0; k < numBins; ++k) {
                    phaseAccum[k] = std::atan2(im[k], re[k]);
                    prevPhase[k]  = phaseAccum[k];
                }
            } else {
                for (size_t k = 0; k < numBins; ++k) {
                    const double mag      = std::sqrt(re[k] * re[k] + im[k] * im[k]);
                    const double phase    = std::atan2(im[k], re[k]);
                    const double expected = k_2pi * static_cast<double>(k) * Ha / N;
                    const double delta    = wrapPhase(phase - prevPhase[k] - expected);

                    phaseAccum[k] += (expected + delta) / static_cast<double>(Ha)
                                     * static_cast<double>(Hs);
                    prevPhase[k]  = phase;

                    re[k] = mag * std::cos(phaseAccum[k]);
                    im[k] = mag * std::sin(phaseAccum[k]);
                }
            }

            // Reconstruct synthesis frame for first frame (phases already set above)
            if (f == 0) {
                for (size_t k = 0; k < numBins; ++k) {
                    const double mag = std::sqrt(re[k] * re[k] + im[k] * im[k]);
                    re[k] = mag * std::cos(phaseAccum[k]);
                    im[k] = mag * std::sin(phaseAccum[k]);
                }
            }

            // Hermitian symmetry → real-valued IFFT output
            for (size_t k = 1; k + 1 < N; ++k) {
                re[N - k] =  re[k];
                im[N - k] = -im[k];
            }

            plan.inverseTransform(re, im, /*scale=*/true);

            // --- Synthesis: windowed overlap-add ---
            for (size_t i = 0; i < N && pos_s + i < outBufLen; ++i)
                outBuf[(pos_s + i) * numChannels + ch] +=
                    static_cast<float>(re[i]) * win[i];
        }
    }

    // Normalise by the window OLA sum and trim to the target length.
    const size_t trimFrames = std::min(static_cast<size_t>(outputFrames), outBufLen);
    for (size_t i = 0; i < trimFrames; ++i) {
        const double norm = normBuf[i] > 1e-10 ? normBuf[i] : 1.0;
        for (int ch = 0; ch < numChannels; ++ch)
            outBuf[i * numChannels + ch] =
                static_cast<float>(outBuf[i * numChannels + ch] / norm);
    }

    Wave result(static_cast<uint32_t>(trimFrames), wave.getNumChannels(),
                wave.getSampleBitDepth(), wave.getSampleRate());
    result.setAudio(outBuf.data(), static_cast<uint32_t>(trimFrames * numChannels));
    wave = std::move(result);
}

} // namespace wm::dsp
