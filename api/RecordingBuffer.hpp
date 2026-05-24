#pragma once

#ifndef RECORDING_BUFFER_H
#define RECORDING_BUFFER_H

#include "Wave.hpp"
#include <atomic>
#include <cstdint>
#include <vector>

namespace wm {

// Accumulates interleaved float samples from a real-time audio input stream.
// append() is called on the audio thread; toWave() is called on the main thread
// after the stream has been stopped (no concurrent access to m_samples).
class RecordingBuffer {
    std::vector<float>    m_samples;
    std::atomic<uint32_t> m_frameCount { 0 };
    uint16_t              m_numChannels = 2;
    uint32_t              m_sampleRate  = 44100;
    uint16_t              m_bitDepth    = 16;

public:
    RecordingBuffer() = default;

    void reset(uint16_t numChannels, uint32_t sampleRate, uint16_t bitDepth,
               uint32_t reserveFrames = 0)
    {
        m_numChannels = numChannels;
        m_sampleRate  = sampleRate;
        m_bitDepth    = bitDepth;
        m_samples.clear();
        if (reserveFrames)
            m_samples.reserve(static_cast<size_t>(reserveFrames) * numChannels);
        m_frameCount.store(0, std::memory_order_relaxed);
    }

    // Audio thread only.
    void append(const float* data, uint32_t numFrames)
    {
        const size_t n = static_cast<size_t>(numFrames) * m_numChannels;
        m_samples.insert(m_samples.end(), data, data + n);
        m_frameCount.fetch_add(numFrames, std::memory_order_release);
    }

    // Safe to read from any thread.
    uint32_t numFrames() const
    {
        return m_frameCount.load(std::memory_order_acquire);
    }

    bool isEmpty() const
    {
        return m_frameCount.load(std::memory_order_relaxed) == 0;
    }

    // Safe to call from the main thread while recording is active.
    // Reads only frames fully committed by the audio thread (up to m_frameCount).
    // Requires the vector to have been pre-reserved so no reallocation occurs.
    Wave peekWave() const
    {
        const uint32_t nf = m_frameCount.load(std::memory_order_acquire);
        if (nf == 0) return Wave();
        Wave wav(nf, m_numChannels, m_bitDepth, m_sampleRate);
        wav.setAudio(m_samples.data(), nf * m_numChannels);
        return wav;
    }

    // Main thread only, call after the recording stream is stopped.
    Wave toWave() const
    {
        if (m_samples.empty())
            return Wave();
        const uint32_t nf = static_cast<uint32_t>(m_samples.size()) / m_numChannels;
        Wave wav(nf, m_numChannels, m_bitDepth, m_sampleRate);
        wav.setAudio(m_samples.data(), static_cast<uint32_t>(m_samples.size()));
        return wav;
    }
};

} // namespace wm

#endif // !RECORDING_BUFFER_H
