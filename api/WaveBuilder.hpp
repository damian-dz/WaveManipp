#pragma once

#ifndef WAVE_BUILDER_H
#define WAVE_BUILDER_H

#include "Wave.hpp"

namespace wm {
 /*!
 * \brief A class that is useful for joining many Wave objects together.
 * \details This class makes it possible to join multiple Wave objects without memory reallocation overhead.
            WaveBuilder objects store pointers to the actual Wave objects.
 */
class WaveBuilder
{
    std::vector<const Wave*> m_wavPtrs;

    uint16_t m_bitsPerSample;
    uint16_t m_numChannels;
    uint32_t m_sampleRate;

public:
    WAVEMANIPPAPI WaveBuilder(uint16_t bitsPerSample, uint16_t numChannels, uint32_t sampleRate);
    WAVEMANIPPAPI WaveBuilder(const Wave& wav);
    WAVEMANIPPAPI WaveBuilder(const WaveBuilder& other);
    WAVEMANIPPAPI ~WaveBuilder();

    WAVEMANIPPAPI WaveBuilder& append(const Wave& wav);
    WAVEMANIPPAPI uint16_t getBitsPerSample() const;
    WAVEMANIPPAPI uint16_t getNumChannels() const;
    WAVEMANIPPAPI uint32_t getSampleRate() const;
    WAVEMANIPPAPI int getWaveCount() const;
    WAVEMANIPPAPI void insert(int idx, const Wave& wav);
    WAVEMANIPPAPI bool isEmpty() const;
    WAVEMANIPPAPI void reserve(uint32_t count);

    WAVEMANIPPAPI Wave toWave();
};
}

#endif // !WAVE_BUILDER_H
