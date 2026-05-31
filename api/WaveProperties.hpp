#pragma once

#ifndef WAVE_PROPERTIES_H
#define WAVE_PROPERTIES_H

#include "utils.hpp"

namespace wm {
/*!
 * \brief Stores RIFF/RIFX WAVE format metadata.
 *
 * \details WaveProperties is a plain metadata container. Setters do not derive or
 * update related values, so callers are responsible for keeping fields such as byte
 * rate, block alignment, data size, and RIFF size consistent. Endianness is tracked
 * by Wave, not by this class.
 */
class WAVEMANIPPAPI WaveProperties
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
    /*! \brief Creates a zero-initialized metadata set. */
    WaveProperties();

    /*! \brief Sets the RIFF chunk size in bytes. */
    void setRiffChunkSize(uint32_t riffChunkSize);
    /*! \brief Sets the fmt chunk size in bytes. */
    void setFmtChunkSize(uint32_t fmtChunkSize);
    /*! \brief Sets the WAVE audio format code, for example 1 for PCM or 3 for float. */
    void setAudioFormat(uint16_t audioFormat);
    /*! \brief Sets the number of interleaved channels. */
    void setNumChannels(uint16_t numChannels);
    /*! \brief Sets the sample rate in Hz. */
    void setSamplingFrequency(uint32_t samplingFrequency);
    /*! \brief Sets the byte rate field. */
    void setNumBytesPerSecond(uint32_t numBytesPerSecond);
    /*! \brief Sets the byte size of one interleaved audio frame. */
    void setBlockAlign(uint16_t blockAlign);
    /*! \brief Sets the number of bits per sample. */
    void setNumBitsPerSample(uint16_t numBitsPerSample);
    /*! \brief Sets the data chunk payload size in bytes. */
    void setDataChunkSize(uint32_t dataChunkSize);

    /*! \brief Returns the RIFF chunk size in bytes. */
    uint32_t getRiffChunkSize() const;
    /*! \brief Returns the fmt chunk size in bytes. */
    uint32_t getFmtChunkSize() const;
    /*! \brief Returns the WAVE audio format code. */
    uint16_t getAudioFormat() const;
    /*! \brief Returns the number of interleaved channels. */
    uint16_t getNumChannels() const;
    /*! \brief Returns the sample rate in Hz. */
    uint32_t getSamplingFrequency() const;
    /*! \brief Returns the byte rate field. */
    uint32_t getNumBytesPerSecond() const;
    /*! \brief Returns the byte size of one interleaved audio frame. */
    uint16_t getBlockAlign() const;
    /*! \brief Returns the number of bits per sample. */
    uint16_t getNumBitsPerSample() const;
    /*! \brief Returns the data chunk payload size in bytes. */
    uint32_t getDataChunkSize() const;
};
}

#endif // !WAVE_PROPERTIES_H
