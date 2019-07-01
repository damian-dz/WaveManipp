#include "WaveProperties.hpp"

namespace wm {

    WaveProperties::WaveProperties() :
        m_riffChunkSize(0),
        m_fmtChunkSize(0),
        m_audioFormat(0),
        m_numChannels(0),
        m_samplingFrequency(0),
        m_numBytesPerSecond(0),
        m_numBitsPerSample(0),
        m_dataChunkSize(0)
    {

    }

    void WaveProperties::setRiffChunkSize(uint32_t riffChunkSize)
    {
        m_riffChunkSize = riffChunkSize;
    }

    void WaveProperties::setFmtChunkSize(uint32_t fmtChunkSize)
    {
        m_fmtChunkSize = fmtChunkSize;
    }

    void WaveProperties::setAudioFormat(uint16_t audioFormat)
    {
        m_audioFormat = audioFormat;
    }

    void WaveProperties::setNumChannels(uint16_t numChannels)
    {
        m_numChannels = numChannels;
    }

    void WaveProperties::setSamplingFrequency(uint32_t samplingFrequency)
    {
        m_samplingFrequency = samplingFrequency;
    }

    void WaveProperties::setNumBytesPerSecond(uint32_t numBytesPerSecond)
    {
        m_numBytesPerSecond = numBytesPerSecond;
    }

    void WaveProperties::setBlockAlign(uint16_t blockAlign)
    {
        m_blockAlign = blockAlign;
    }

    void WaveProperties::setNumBitsPerSample(uint16_t numBitsPerSample)
    {
        m_numBitsPerSample = numBitsPerSample;
    }

    void WaveProperties::setDataChunkSize(uint32_t dataChunkSize)
    {
        m_dataChunkSize = dataChunkSize;
    }

    uint32_t WaveProperties::getRiffChunkSize()
    {
        return m_riffChunkSize;
    }

    uint32_t WaveProperties::getFmtChunkSize()
    {
        return m_fmtChunkSize;
    }

    uint16_t WaveProperties::getAudioFormat()
    {
        return m_audioFormat;
    }

    uint16_t WaveProperties::getNumChannels()
    {
        return m_numChannels;
    }

    uint32_t WaveProperties::getSamplingFrequency()
    {
        return m_samplingFrequency;
    }

    uint32_t WaveProperties::getNumBytesPerSecond()
    {
        return m_numBytesPerSecond;
    }

    uint16_t WaveProperties::getBlockAlign()
    {
        return m_blockAlign;
    }

    uint16_t WaveProperties::getNumBitsPerSample()
    {
        return m_numBitsPerSample;
    }

    uint32_t WaveProperties::getDataChunkSize()
    {
        return m_dataChunkSize;
    }

}
