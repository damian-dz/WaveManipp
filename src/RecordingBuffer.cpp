#include "RecordingBuffer.hpp"

namespace wm {

void RecordingBuffer::reset(uint16_t numChannels, uint32_t sampleRate, uint16_t bitDepth,
                            uint32_t reserveFrames)
{
    m_numChannels = numChannels;
    m_sampleRate  = sampleRate;
    m_bitDepth    = bitDepth;
    m_samples.clear();
    if (reserveFrames)
        m_samples.reserve(static_cast<size_t>(reserveFrames) * numChannels);
    m_frameCount.store(0, std::memory_order_relaxed);
}

void RecordingBuffer::append(const float* data, uint32_t numFrames)
{
    const size_t n = static_cast<size_t>(numFrames) * m_numChannels;
    m_samples.insert(m_samples.end(), data, data + n);
    m_frameCount.fetch_add(numFrames, std::memory_order_release);
}

uint32_t RecordingBuffer::numFrames() const
{
    return m_frameCount.load(std::memory_order_acquire);
}

bool RecordingBuffer::isEmpty() const
{
    return m_frameCount.load(std::memory_order_relaxed) == 0;
}

Wave RecordingBuffer::peekWave() const
{
    const uint32_t nf = m_frameCount.load(std::memory_order_acquire);
    if (nf == 0) return Wave();
    Wave wav(nf, m_numChannels, m_bitDepth, m_sampleRate);
    wav.setAudio(m_samples.data(), nf * m_numChannels);
    return wav;
}

Wave RecordingBuffer::toWave() const
{
    if (m_samples.empty())
        return Wave();
    const uint32_t nf = static_cast<uint32_t>(m_samples.size()) / m_numChannels;
    Wave wav(nf, m_numChannels, m_bitDepth, m_sampleRate);
    wav.setAudio(m_samples.data(), static_cast<uint32_t>(m_samples.size()));
    return wav;
}

} // namespace wm
