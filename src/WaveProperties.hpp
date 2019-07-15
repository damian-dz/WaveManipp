#pragma once

#ifndef WAVE_PROPERTIES_H
#define WAVE_PROPERTIES_H

#include "utils.hpp"

namespace wm {
/*!
 * \brief A class that encapsulates the typical properties of a RIFF/RIFX WAVE file.
 * \details None of the setter methods in this class performs any sort of arithmetic,
            so setting one of the parameters will not affect any other.
            The class stores no information about endianness.
 */
class WaveProperties
{
    uint32_t m_riffChunkSize;
    uint32_t m_fmtChunkSize;
    uint16_t m_audioFormat;
    uint16_t m_numChannels;
    uint32_t m_samplingFrequency;
    uint32_t m_numBytesPerSecond;
    uint16_t m_blockAlign;
    uint16_t m_numBitsPerSample;
    uint32_t m_dataChunkSize;

public:
    WaveProperties();

    void setRiffChunkSize(uint32_t riffChunkSize);
    void setFmtChunkSize(uint32_t fmtChunkSize);
    void setAudioFormat(uint16_t audioFormat);
    void setNumChannels(uint16_t numChannels);
    void setSamplingFrequency(uint32_t samplingFrequency);
    void setNumBytesPerSecond(uint32_t numBytesPerSecond);
    void setBlockAlign(uint16_t blockAlign);
    void setNumBitsPerSample(uint16_t numBitsPerSample);
    void setDataChunkSize(uint32_t dataChunkSize);

    uint32_t getRiffChunkSize() const;
    uint32_t getFmtChunkSize() const;
    uint16_t getAudioFormat() const;
    uint16_t getNumChannels() const;
    uint32_t getSamplingFrequency() const;
    uint32_t getNumBytesPerSecond() const;
    uint16_t getBlockAlign() const;
    uint16_t getNumBitsPerSample() const;
    uint32_t getDataChunkSize() const;
};
}

#endif // !WAVE_PROPERTIES_H
