#include "WaveProperties.hpp"

namespace wm {
/*!
 * \brief Creates an empty (zero-initialized) instance of WaveProperties.
 */
WaveProperties::WaveProperties() :
    m_riffChunkSize(0),
    m_fmtChunkSize(0),
    m_audioFormat(0),
    m_numChannels(0),
    m_samplingFrequency(0),
    m_numBytesPerSecond(0),
    m_blockAlign(0),
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

/*!
 * \brief Sets the sampling frequency for the audio data.
 * \details Also known as <i>sampling rate</i> and <i>sample rate</i>.
 * \param samplingFrequency &mdash; the new sampling frequency
 */
void WaveProperties::setSamplingFrequency(uint32_t samplingFrequency)
{
    m_samplingFrequency = samplingFrequency;
}

/*!
 * \brief Sets the number of bytes necessary to store one second of the audio data.
 * \details Also known as <i>byte rate</i>.
 * \param numBytesPerSecond &mdash; the new number of bytes per second
 */
void WaveProperties::setNumBytesPerSecond(uint32_t numBytesPerSecond)
{
    m_numBytesPerSecond = numBytesPerSecond;
}

void WaveProperties::setBlockAlign(uint16_t blockAlign)
{
    m_blockAlign = blockAlign;
}

/*!
 * \brief Sets the number of bits required to represent each audio sample.
 * \details Also known as <i>sample bit depth</i>.
 * \param numBitsPerSample &mdash; the new number of bits per sample
 */
void WaveProperties::setNumBitsPerSample(uint16_t numBitsPerSample)
{
    m_numBitsPerSample = numBitsPerSample;
}

void WaveProperties::setDataChunkSize(uint32_t dataChunkSize)
{
    m_dataChunkSize = dataChunkSize;
}

uint32_t WaveProperties::getRiffChunkSize() const
{
    return m_riffChunkSize;
}

uint32_t WaveProperties::getFmtChunkSize() const
{
    return m_fmtChunkSize;
}

uint16_t WaveProperties::getAudioFormat() const
{
    return m_audioFormat;
}

uint16_t WaveProperties::getNumChannels() const
{
    return m_numChannels;
}

uint32_t WaveProperties::getSamplingFrequency() const
{
    return m_samplingFrequency;
}

uint32_t WaveProperties::getNumBytesPerSecond() const
{
    return m_numBytesPerSecond;
}

uint16_t WaveProperties::getBlockAlign() const
{
    return m_blockAlign;
}

uint16_t WaveProperties::getNumBitsPerSample() const
{
    return m_numBitsPerSample;
}

uint32_t WaveProperties::getDataChunkSize() const
{
    return m_dataChunkSize;
}

}
