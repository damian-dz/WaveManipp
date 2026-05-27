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

} // namespace wm::dsp
