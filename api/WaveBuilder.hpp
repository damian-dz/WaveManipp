#pragma once

/*!
 * \file WaveBuilder.hpp
 * \brief WaveBuilder class for sequentially concatenating compatible Wave objects.
 */

#ifndef WAVE_BUILDER_H
#define WAVE_BUILDER_H

#include "Wave.hpp"

namespace wm {
 /*!
 * \brief Concatenates multiple compatible Wave objects into one Wave.
 *
 * \details WaveBuilder stores non-owning pointers to the source Wave objects and
 * copies their audio only when toWave() is called. All appended waves must have the
 * same bit depth, channel count, and sample rate. The source Wave instances must
 * outlive the builder or at least remain valid until toWave() finishes.
 */
class WaveBuilder
{
    std::vector<const Wave*> m_wavPtrs;

    uint16_t m_bitsPerSample;
    uint16_t m_numChannels;
    uint32_t m_sampleRate;

public:
    /*!
     * \brief Creates an empty builder for waves with the requested format.
     * \param bitsPerSample Required bit depth of all appended waves.
     * \param numChannels   Required channel count of all appended waves.
     * \param sampleRate    Required sample rate in Hz of all appended waves.
     */
    WAVEMANIPPAPI WaveBuilder(uint16_t bitsPerSample, uint16_t numChannels, uint32_t sampleRate);
    /*! \brief Creates a builder initialized with one source wave. */
    WAVEMANIPPAPI WaveBuilder(const Wave& wav);
    /*! \brief Copies the list of non-owning Wave pointers and format settings. */
    WAVEMANIPPAPI WaveBuilder(const WaveBuilder& other);
    /*! \brief Clears the builder's pointer list. */
    WAVEMANIPPAPI ~WaveBuilder();

    /*!
     * \brief Appends a compatible Wave to the end of the output sequence.
     * \param wav Source wave; must match this builder's bit depth, channel count, and sample rate.
     * \return Reference to this builder for method chaining.
     */
    WAVEMANIPPAPI WaveBuilder& append(const Wave& wav);
    /*! \brief Returns the required source bit depth. */
    WAVEMANIPPAPI uint16_t getBitsPerSample() const;
    /*! \brief Returns the required source channel count. */
    WAVEMANIPPAPI uint16_t getNumChannels() const;
    /*! \brief Returns the required source sample rate in Hz. */
    WAVEMANIPPAPI uint32_t getSampleRate() const;
    /*! \brief Returns the number of source waves currently referenced. */
    WAVEMANIPPAPI int getWaveCount() const;
    /*!
     * \brief Inserts a compatible Wave before the element at the specified index.
     * \param idx Zero-based insertion position; existing elements shift right.
     * \param wav Source wave; must match this builder's bit depth, channel count, and sample rate.
     */
    WAVEMANIPPAPI void insert(int idx, const Wave& wav);
    /*! \brief Returns true when no source waves are referenced. */
    WAVEMANIPPAPI bool isEmpty() const;
    /*!
     * \brief Reserves pointer-list capacity for a number of source waves.
     * \param count Number of wave pointers to pre-allocate capacity for.
     */
    WAVEMANIPPAPI void reserve(uint32_t count);

    /*! \brief Copies all referenced waves into a single concatenated Wave. */
    WAVEMANIPPAPI Wave toWave();
};
}

#endif // !WAVE_BUILDER_H
