#pragma once

/*!
 * \file Dsp.hpp
 * \brief In-place DSP effects, analysis, and filtering for Wave objects.
 */

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
/*!
 * \brief Returns the largest absolute sample value within a frame range.
 * \param wave       Source wave to analyze.
 * \param startFrame First frame of the range (inclusive).
 * \param endFrame   One past the last frame of the range (exclusive).
 */
float peakRange(const Wave& wave, uint32_t startFrame, uint32_t endFrame);

// ---- Gain -------------------------------------------------------------------

/*!
 * \brief Multiplies every sample by a linear gain value.
 * \param wave        Wave to process in place.
 * \param linearGain  Linear gain factor (1.0 = unity).
 */
void applyGain(Wave& wave, float linearGain);
/*!
 * \brief Multiplies every sample in a frame range by a linear gain value.
 * \param wave        Wave to process in place.
 * \param startFrame  First frame of the range (inclusive).
 * \param endFrame    One past the last frame of the range (exclusive).
 * \param linearGain  Linear gain factor (1.0 = unity).
 */
void applyGainRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float linearGain);

/*!
 * \brief Normalizes the whole wave to a target peak level in dBFS.
 * \param wave     Wave to normalize in place.
 * \param targetDb Target peak level in dBFS (0.0 = full scale).
 */
void normalize(Wave& wave, float targetDb = 0.f);
/*!
 * \brief Normalizes a frame range to a target peak level in dBFS.
 * \param wave       Wave to normalize in place.
 * \param startFrame First frame of the range (inclusive).
 * \param endFrame   One past the last frame of the range (exclusive).
 * \param targetDb   Target peak level in dBFS (0.0 = full scale).
 */
void normalizeRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, float targetDb = 0.f);

// ---- Editing ----------------------------------------------------------------

/*!
 * \brief Sets every sample in a frame range to zero.
 * \param wave       Wave to silence in place.
 * \param startFrame First frame of the range (inclusive).
 * \param endFrame   One past the last frame of the range (exclusive).
 */
void silence(Wave& wave, uint32_t startFrame, uint32_t endFrame);

/*! \brief Changes the playback speed by \p factor (pitch shifts proportionally).
 *
 *  The wave is reinterpreted at \p factor × its current sample rate, then
 *  resampled back to the original rate. Factor > 1 makes it shorter and
 *  higher-pitched; factor < 1 makes it longer and lower-pitched.
 *  \param wave   Wave to process in place.
 *  \param factor Speed multiplier (> 1 = shorter/higher, < 1 = longer/lower).
 */
void changeSpeed(Wave& wave, float factor);

/*! \brief Inverts polarity for the whole wave. \param wave Wave to invert in place. */
void invert(Wave& wave);
/*!
 * \brief Inverts polarity for a frame range.
 * \param wave       Wave to invert in place.
 * \param startFrame First frame of the range (inclusive).
 * \param endFrame   One past the last frame of the range (exclusive).
 */
void invertRange(Wave& wave, uint32_t startFrame, uint32_t endFrame);

// ---- Fades ------------------------------------------------------------------

/*! \brief Fade curve shape used by fadeIn() and fadeOut(). */
enum class FadeCurve { Linear, Logarithmic, EqualPower, SCurve };

/*!
 * \brief Applies a fade-in over a frame range.
 * \param wave       Wave to process in place.
 * \param startFrame First frame of the fade (inclusive).
 * \param endFrame   One past the last frame of the fade (exclusive).
 * \param curve      Shape of the fade envelope.
 */
void fadeIn (Wave& wave, uint32_t startFrame, uint32_t endFrame, FadeCurve curve = FadeCurve::Logarithmic);
/*!
 * \brief Applies a fade-out over a frame range.
 * \param wave       Wave to process in place.
 * \param startFrame First frame of the fade (inclusive).
 * \param endFrame   One past the last frame of the fade (exclusive).
 * \param curve      Shape of the fade envelope.
 */
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

/*!
 * \brief Builds a window coefficient table.
 * \param windowSize Number of coefficients to generate.
 * \param window     Window function type.
 * \return Vector of \p windowSize normalized window coefficients.
 */
std::vector<float> makeWindow(size_t windowSize, WindowFunction window);
/*!
 * \brief Computes single-sided STFT magnitude frames in dB.
 * \param wave   Wave to analyze.
 * \param config Analysis parameters (window size, hop size, channel, floor).
 * \return 2-D array of magnitude values in dB: [frame][bin].
 */
std::vector<std::vector<float>> stftMagnitudeDb(const Wave& wave, const StftConfig& config = {});

// ---- Reverb -----------------------------------------------------------------

/*! \brief Parameters for the Schroeder/Freeverb-style reverb. */
struct ReverbParams {
    float roomSize = 0.5f;  ///< Comb feedback amount (0 = tiny room, 1 = large hall)
    float damping  = 0.5f;  ///< High-frequency absorption (0 = bright, 1 = dark)
    float wetMix   = 0.3f;  ///< Dry/wet blend (0 = all dry, 1 = all wet)
};

/*!
 * \brief Applies reverb over a frame range. The tail is clipped at endFrame.
 * \param wave       Wave to process in place.
 * \param startFrame First frame of the region to reverberate (inclusive).
 * \param endFrame   One past the last frame; reverb tail is clipped here (exclusive).
 * \param params     Room size, damping, and wet-mix parameters.
 */
void reverb(Wave& wave, uint32_t startFrame, uint32_t endFrame, const ReverbParams& params = {});

// ---- Tempo / Pitch ----------------------------------------------------------

/*! \brief Stretches or compresses the wave in time without changing pitch.
 *
 *  Uses a phase vocoder (Hann-windowed STFT with phase accumulation).
 *  \p factor > 1 makes the wave longer (slower tempo); < 1 makes it shorter (faster tempo).
 *  \p fftSize must be a power of two; 2048 is a good default.
 *  \param wave    Wave to stretch in place.
 *  \param factor  Time-stretch factor (> 1 = slower/longer, < 1 = faster/shorter).
 *  \param fftSize FFT window size in samples; must be a power of two.
 */
void stretchTempo(Wave& wave, float factor, size_t fftSize = 2048);

// ---- Reverse ----------------------------------------------------------------

/*!
 * \brief Reverses one channel in place.
 * \param wave    Wave to process in place.
 * \param channel Zero-based channel index to reverse.
 */
void reverseChannel(Wave& wave, int channel);
/*! \brief Reverses every channel in place. \param wave Wave to reverse in place. */
void reverse(Wave& wave);

// ---- Biquad IIR filter ------------------------------------------------------

/*! \brief Normalized Direct Form I biquad coefficients with a0 assumed to be 1. */
struct BiquadCoeffs {
    float b0, b1, b2;   ///< Feed-forward coefficients (normalized, a0 = 1)
    float a1, a2;        ///< Feed-back coefficients    (normalized, a0 = 1)
};

/*!
 * \brief Creates coefficients for a peaking (bell) EQ filter.
 * \param freqHz     Center frequency in Hz.
 * \param sampleRate Sample rate of the audio to be processed, in Hz.
 * \param gainDb     Boost or cut in dB (positive = boost, negative = cut).
 * \param Q          Quality factor controlling bandwidth (higher Q = narrower).
 */
BiquadCoeffs makePeaking  (float freqHz, float sampleRate, float gainDb, float Q);
/*!
 * \brief Creates coefficients for a low-shelf EQ filter.
 * \param freqHz     Shelf transition frequency in Hz.
 * \param sampleRate Sample rate of the audio to be processed, in Hz.
 * \param gainDb     Boost or cut in dB for frequencies below \p freqHz.
 */
BiquadCoeffs makeLowShelf (float freqHz, float sampleRate, float gainDb);
/*!
 * \brief Creates coefficients for a high-shelf EQ filter.
 * \param freqHz     Shelf transition frequency in Hz.
 * \param sampleRate Sample rate of the audio to be processed, in Hz.
 * \param gainDb     Boost or cut in dB for frequencies above \p freqHz.
 */
BiquadCoeffs makeHighShelf(float freqHz, float sampleRate, float gainDb);

/*!
 * \brief Processes the whole wave with a biquad filter.
 * \param wave Wave to filter in place.
 * \param c    Normalized biquad coefficients.
 */
void processBiquad     (Wave& wave, const BiquadCoeffs& c);
/*!
 * \brief Processes a frame range with a biquad filter.
 * \param wave       Wave to filter in place.
 * \param startFrame First frame of the range (inclusive).
 * \param endFrame   One past the last frame of the range (exclusive).
 * \param c          Normalized biquad coefficients.
 */
void processBiquadRange(Wave& wave, uint32_t startFrame, uint32_t endFrame, const BiquadCoeffs& c);

} // namespace wm::dsp
