#pragma once
#include "Wave.hpp"

namespace wm::dsp {

// ---- Queries ----------------------------------------------------------------

float peak(const Wave& wave);
float peakRange(const Wave& wave, uint32_t startFrame, uint32_t endFrame);

// ---- Gain -------------------------------------------------------------------

void applyGain(Wave& wave, float linearGain);
void applyGainRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float linearGain);

void normalize(Wave& wave, float targetDb = 0.f);
void normalizeRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float targetDb = 0.f);

// ---- Editing ----------------------------------------------------------------

void silence(Wave& wave, uint32_t startFrame, uint32_t endFrame);

void invert(Wave& wave);
void invertRange(Wave& wave, uint32_t startFrame, uint32_t endFrame);

// ---- Fades ------------------------------------------------------------------

enum class FadeCurve { Linear, Logarithmic };

void fadeIn (Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve = FadeCurve::Logarithmic);
void fadeOut(Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve = FadeCurve::Logarithmic);

// ---- Short-time Fourier transform ------------------------------------------

enum class WindowFunction {
    Rectangular,
    Hann,
    Hamming,
    Blackman
};

struct StftConfig {
    size_t windowSize = 2048;
    size_t hopSize = 512;
    WindowFunction window = WindowFunction::Hann;
    int channel = 0;
    float floorDb = -120.f;
};

std::vector<float> makeWindow(size_t windowSize, WindowFunction window);
std::vector<std::vector<float>> stftMagnitudeDb(const Wave& wave, const StftConfig& config = {});

// ---- Reverse ----------------------------------------------------------------

void reverseChannel(Wave& wave, int channel);
void reverse(Wave& wave);

// ---- Biquad IIR filter ------------------------------------------------------

struct BiquadCoeffs {
    float b0, b1, b2;   // feed-forward coefficients (normalized, a0 = 1)
    float a1, a2;        // feed-back coefficients    (normalized, a0 = 1)
};

BiquadCoeffs makePeaking  (float freqHz, float sampleRate, float gainDb, float Q);
BiquadCoeffs makeLowShelf (float freqHz, float sampleRate, float gainDb);
BiquadCoeffs makeHighShelf(float freqHz, float sampleRate, float gainDb);

void processBiquad     (Wave& wave, const BiquadCoeffs& c);
void processBiquadRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, const BiquadCoeffs& c);

} // namespace wm::dsp
