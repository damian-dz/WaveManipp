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
    WaveBuilder(uint16_t bitsPerSample, uint16_t numChannels, uint32_t sampleRate);
    WaveBuilder(const Wave& wav);
    WaveBuilder(const WaveBuilder& other);
    ~WaveBuilder();

    WaveBuilder& append(const Wave& wav);
    uint16_t getBitsPerSample() const;
    uint16_t getNumChannels() const;
    uint32_t getSampleRate() const;
    int getWaveCount() const;
    void insert(int idx, const Wave& wav);
    bool isEmpty() const;
    void reserve(uint32_t count);

    Wave toWave();
};
}

#endif // !WAVE_BUILDER_H
