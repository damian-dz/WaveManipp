#include "Dsp.hpp"
#include <cmath>
#include <algorithm>

namespace wm::dsp {

float peak(const Wave& wave)
{
    if (wave.isEmpty()) return 0.f;
    const uint16_t nc = wave.getNumChannels();
    const uint32_t nf = wave.getNumFrames();
    float p = 0.f;
    for (uint32_t fr = 0; fr < nf; ++fr)
        for (uint16_t ch = 0; ch < nc; ++ch)
            p = std::max(p, std::fabs(wave(fr, static_cast<int>(ch))));
    return p;
}

float peakRange(const Wave& wave, uint32_t startFrame, uint32_t endFrame)
{
    if (wave.isEmpty() || startFrame >= endFrame) return 0.f;
    const uint16_t nc  = wave.getNumChannels();
    const uint32_t end = std::min(endFrame, wave.getNumFrames());
    float p = 0.f;
    for (uint32_t fr = startFrame; fr < end; ++fr)
        for (uint16_t ch = 0; ch < nc; ++ch)
            p = std::max(p, std::fabs(wave(fr, static_cast<int>(ch))));
    return p;
}

void applyGain(Wave& wave, float linearGain)
{
    const uint16_t nc = wave.getNumChannels();
    const uint32_t nf = wave.getNumFrames();
    for (uint32_t fr = 0; fr < nf; ++fr)
        for (uint16_t ch = 0; ch < nc; ++ch)
            wave(fr, static_cast<int>(ch)) *= linearGain;
}

void applyGainRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float linearGain)
{
    const uint16_t nc  = wave.getNumChannels();
    const uint32_t end = std::min(endFrame, wave.getNumFrames());
    for (uint32_t fr = startFrame; fr < end; ++fr)
        for (uint16_t ch = 0; ch < nc; ++ch)
            wave(fr, static_cast<int>(ch)) *= linearGain;
}

void normalize(Wave& wave, float targetDb)
{
    const float p = peak(wave);
    if (p == 0.f) return;
    applyGain(wave, std::pow(10.f, targetDb / 20.f) / p);
}

void normalizeRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float targetDb)
{
    const float p = peakRange(wave, startFrame, endFrame);
    if (p == 0.f) return;
    applyGainRange(wave, startFrame, endFrame, std::pow(10.f, targetDb / 20.f) / p);
}

void silence(Wave& wave, uint32_t startFrame, uint32_t endFrame)
{
    const uint16_t nc  = wave.getNumChannels();
    const uint32_t end = std::min(endFrame, wave.getNumFrames());
    for (uint32_t fr = startFrame; fr < end; ++fr)
        for (uint16_t ch = 0; ch < nc; ++ch)
            wave(fr, static_cast<int>(ch)) = 0.f;
}

void invert(Wave& wave)
{
    applyGain(wave, -1.f);
}

void invertRange(Wave& wave, uint32_t startFrame, uint32_t endFrame)
{
    applyGainRange(wave, startFrame, endFrame, -1.f);
}

void fadeIn(Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve)
{
    const uint16_t nc  = wave.getNumChannels();
    const uint32_t end = std::min(endFrame, wave.getNumFrames());
    if (startFrame >= end) return;
    const float len = static_cast<float>(end - startFrame);
    for (uint32_t fr = startFrame; fr < end; ++fr) {
        const float t    = static_cast<float>(fr - startFrame) / len;
        const float gain = (curve == FadeCurve::Logarithmic)
            ? std::pow(10.f, -3.f * (1.f - t))  // 0 dB at t=1, -60 dB at t=0
            : t;
        for (uint16_t ch = 0; ch < nc; ++ch)
            wave(fr, static_cast<int>(ch)) *= gain;
    }
}

void fadeOut(Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve)
{
    const uint16_t nc  = wave.getNumChannels();
    const uint32_t end = std::min(endFrame, wave.getNumFrames());
    if (startFrame >= end) return;
    const float len = static_cast<float>(end - startFrame);
    for (uint32_t fr = startFrame; fr < end; ++fr) {
        const float t    = static_cast<float>(fr - startFrame) / len;
        const float gain = (curve == FadeCurve::Logarithmic)
            ? std::pow(10.f, -3.f * t)           // 0 dB at t=0, -60 dB at t=1
            : (1.f - t);
        for (uint16_t ch = 0; ch < nc; ++ch)
            wave(fr, static_cast<int>(ch)) *= gain;
    }
}

void reverseChannel(Wave& wave, int channel)
{
    int lo = 0;
    int hi = static_cast<int>(wave.getNumFrames()) - 1;
    while (lo < hi) {
        std::swap(wave(static_cast<uint32_t>(lo), channel),
                  wave(static_cast<uint32_t>(hi), channel));
        ++lo; --hi;
    }
}

void reverse(Wave& wave)
{
    const uint16_t nc = wave.getNumChannels();
    for (uint16_t ch = 0; ch < nc; ++ch)
        reverseChannel(wave, static_cast<int>(ch));
}

// ---- Biquad IIR filter ------------------------------------------------------
// Coefficient formulas from "Audio EQ Cookbook" by Robert Bristow-Johnson.
// All filters use S = 1 (maximally flat shelf slope) and Direct Form I.

static float clampFreq(float freqHz, float sampleRate)
{
    return std::max(1.f, std::min(freqHz, sampleRate * 0.5f - 1.f));
}

BiquadCoeffs makePeaking(float freqHz, float sampleRate, float gainDb, float Q)
{
    const float A     = std::pow(10.f, gainDb / 40.f);
    const float w0    = 2.f * 3.14159265358979f * clampFreq(freqHz, sampleRate) / sampleRate;
    const float alpha = std::sin(w0) / (2.f * std::max(Q, 0.001f));
    const float a0    = 1.f + alpha / A;
    return {
        (1.f + alpha * A) / a0,
        (-2.f * std::cos(w0)) / a0,
        (1.f - alpha * A) / a0,
        (-2.f * std::cos(w0)) / a0,
        (1.f - alpha / A) / a0
    };
}

BiquadCoeffs makeLowShelf(float freqHz, float sampleRate, float gainDb)
{
    const float A     = std::pow(10.f, gainDb / 40.f);
    const float w0    = 2.f * 3.14159265358979f * clampFreq(freqHz, sampleRate) / sampleRate;
    const float alpha = std::sin(w0) / std::sqrt(2.f);   // S = 1
    const float sqA   = std::sqrt(A);
    const float cosw  = std::cos(w0);
    const float a0    = (A + 1.f) + (A - 1.f) * cosw + 2.f * sqA * alpha;
    return {
         A * ((A + 1.f) - (A - 1.f) * cosw + 2.f * sqA * alpha) / a0,
         2.f * A * ((A - 1.f) - (A + 1.f) * cosw) / a0,
         A * ((A + 1.f) - (A - 1.f) * cosw - 2.f * sqA * alpha) / a0,
        -2.f * ((A - 1.f) + (A + 1.f) * cosw) / a0,
         ((A + 1.f) + (A - 1.f) * cosw - 2.f * sqA * alpha) / a0
    };
}

BiquadCoeffs makeHighShelf(float freqHz, float sampleRate, float gainDb)
{
    const float A     = std::pow(10.f, gainDb / 40.f);
    const float w0    = 2.f * 3.14159265358979f * clampFreq(freqHz, sampleRate) / sampleRate;
    const float alpha = std::sin(w0) / std::sqrt(2.f);   // S = 1
    const float sqA   = std::sqrt(A);
    const float cosw  = std::cos(w0);
    const float a0    = (A + 1.f) - (A - 1.f) * cosw + 2.f * sqA * alpha;
    return {
         A * ((A + 1.f) + (A - 1.f) * cosw + 2.f * sqA * alpha) / a0,
        -2.f * A * ((A - 1.f) + (A + 1.f) * cosw) / a0,
         A * ((A + 1.f) + (A - 1.f) * cosw - 2.f * sqA * alpha) / a0,
         2.f * ((A - 1.f) - (A + 1.f) * cosw) / a0,
         ((A + 1.f) - (A - 1.f) * cosw - 2.f * sqA * alpha) / a0
    };
}

void processBiquadRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, const BiquadCoeffs& c)
{
    const uint16_t nc  = wave.getNumChannels();
    const uint32_t end = std::min(endFrame, wave.getNumFrames());
    if (startFrame >= end) return;
    for (uint16_t ch = 0; ch < nc; ++ch) {
        float x1 = 0.f, x2 = 0.f, y1 = 0.f, y2 = 0.f;
        for (uint32_t fr = startFrame; fr < end; ++fr) {
            const float x0 = wave(fr, static_cast<int>(ch));
            const float y0 = c.b0*x0 + c.b1*x1 + c.b2*x2 - c.a1*y1 - c.a2*y2;
            wave(fr, static_cast<int>(ch)) = y0;
            x2 = x1; x1 = x0;
            y2 = y1; y1 = y0;
        }
    }
}

void processBiquad(Wave& wave, const BiquadCoeffs& c)
{
    processBiquadRange(wave, 0, wave.getNumFrames(), c);
}

} // namespace wm::dsp
