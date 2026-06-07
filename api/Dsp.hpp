#pragma once
#include "Wave.hpp"

/*!
 * \namespace wm::dsp
 * \brief Utility DSP functions for in-place Wave editing and spectral analysis.
 *
 * \details Frame ranges use half-open intervals: startFrame is included and
 * endFrame is excluded. Range operations clamp endFrame to the Wave length.
 */
namespace wm::dsp {

// ---- Queries ----------------------------------------------------------------

/*! \brief Returns the largest absolute sample value across all channels. */
float peak(const Wave& wave);
/*! \brief Returns the largest absolute sample value in a frame range. */
float peakRange(const Wave& wave, uint32_t startFrame, uint32_t endFrame);

// ---- Gain -------------------------------------------------------------------

/*! \brief Multiplies every sample by a linear gain value. */
void applyGain(Wave& wave, float linearGain);
/*! \brief Multiplies every sample in a frame range by a linear gain value. */
void applyGainRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float linearGain);

/*! \brief Normalizes the whole wave to a target peak level in dBFS. */
void normalize(Wave& wave, float targetDb = 0.f);
/*! \brief Normalizes a frame range to a target peak level in dBFS. */
void normalizeRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float targetDb = 0.f);

// ---- Editing ----------------------------------------------------------------

/*! \brief Sets every sample in a frame range to zero. */
void silence(Wave& wave, uint32_t startFrame, uint32_t endFrame);

/*! \brief Changes the playback speed by \p factor (pitch shifts proportionally).
 *
 *  The wave is reinterpreted at \p factor × its current sample rate, then
 *  resampled back to the original rate. Factor > 1 makes it shorter and
 *  higher-pitched; factor < 1 makes it longer and lower-pitched.
 */
void changeSpeed(Wave& wave, float factor);

/*! \brief Inverts polarity for the whole wave. */
void invert(Wave& wave);
/*! \brief Inverts polarity for a frame range. */
void invertRange(Wave& wave, uint32_t startFrame, uint32_t endFrame);

// ---- Fades ------------------------------------------------------------------

/*! \brief Fade curve shape used by fadeIn() and fadeOut(). */
enum class FadeCurve { Linear, Logarithmic, EqualPower, SCurve };

/*! \brief Applies a fade-in over a frame range. */
void fadeIn (Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve = FadeCurve::Logarithmic);
/*! \brief Applies a fade-out over a frame range. */
void fadeOut(Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve = FadeCurve::Logarithmic);

// ---- Short-time Fourier transform ------------------------------------------

/*! \brief Window functions available for STFT analysis. */
enum class WindowFunction {
    Rectangular,
    Hann,
    Hamming,
    Blackman
};

/*! \brief Configuration for stftMagnitudeDb(). */
struct StftConfig {
    /*! \brief Number of samples per analysis window. */
    size_t windowSize = 2048;
    /*! \brief Number of frames between adjacent analysis windows. */
    size_t hopSize = 512;
    /*! \brief Window shape applied before each FFT. */
    WindowFunction window = WindowFunction::Hann;
    /*! \brief Channel index to analyze. */
    int channel = 0;
    /*! \brief Minimum output level in dB. */
    float floorDb = -120.f;
};

/*! \brief Builds a window coefficient table. */
std::vector<float> makeWindow(size_t windowSize, WindowFunction window);
/*! \brief Computes single-sided STFT magnitude frames in dB. */
std::vector<std::vector<float>> stftMagnitudeDb(const Wave& wave, const StftConfig& config = {});

// ---- Reverb -----------------------------------------------------------------

/*! \brief Parameters for the Schroeder/Freeverb-style reverb. */
struct ReverbParams {
    float roomSize = 0.5f;  ///< Comb feedback amount (0 = tiny room, 1 = large hall)
    float damping  = 0.5f;  ///< High-frequency absorption (0 = bright, 1 = dark)
    float wetMix   = 0.3f;  ///< Dry/wet blend (0 = all dry, 1 = all wet)
};

/*! \brief Applies reverb over a frame range. The tail is clipped at endFrame. */
void reverb(Wave& wave, uint32_t startFrame, uint32_t endFrame, const ReverbParams& params = {});

// ---- Reverse ----------------------------------------------------------------

/*! \brief Reverses one channel in place. */
void reverseChannel(Wave& wave, int channel);
/*! \brief Reverses every channel in place. */
void reverse(Wave& wave);

// ---- Biquad IIR filter ------------------------------------------------------

/*! \brief Normalized Direct Form I biquad coefficients with a0 assumed to be 1. */
struct BiquadCoeffs {
    float b0, b1, b2;   // feed-forward coefficients (normalized, a0 = 1)
    float a1, a2;        // feed-back coefficients    (normalized, a0 = 1)
};

/*! \brief Creates coefficients for a peaking EQ filter. */
BiquadCoeffs makePeaking  (float freqHz, float sampleRate, float gainDb, float Q);
/*! \brief Creates coefficients for a low-shelf EQ filter. */
BiquadCoeffs makeLowShelf (float freqHz, float sampleRate, float gainDb);
/*! \brief Creates coefficients for a high-shelf EQ filter. */
BiquadCoeffs makeHighShelf(float freqHz, float sampleRate, float gainDb);

/*! \brief Processes the whole wave with a biquad filter. */
void processBiquad     (Wave& wave, const BiquadCoeffs& c);
/*! \brief Processes a frame range with a biquad filter. */
void processBiquadRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, const BiquadCoeffs& c);

} // namespace wm::dsp
