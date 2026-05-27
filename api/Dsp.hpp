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

// ---- Reverse ----------------------------------------------------------------

void reverseChannel(Wave& wave, int channel);
void reverse(Wave& wave);

} // namespace wm::dsp
