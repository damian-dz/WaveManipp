#pragma once

#ifndef RECORDING_BUFFER_H
#define RECORDING_BUFFER_H

#include "Wave.hpp"
#include <atomic>
#include <cstdint>
#include <vector>

namespace wm {

// append() is audio-thread only; peekWave() is safe from any thread while recording;
// toWave() is main-thread only after the stream is stopped.
class RecordingBuffer {
    std::vector<float>    m_samples;
    std::atomic<uint32_t> m_frameCount { 0 };
    uint16_t              m_numChannels = 2;
    uint32_t              m_sampleRate  = 44100;
    uint16_t              m_bitDepth    = 16;

public:
    RecordingBuffer() = default;

    WAVEMANIPPAPI void     reset(uint16_t numChannels, uint32_t sampleRate, uint16_t bitDepth,
                                 uint32_t reserveFrames = 0);
    WAVEMANIPPAPI void     append(const float* data, uint32_t numFrames);
    WAVEMANIPPAPI uint32_t numFrames() const;
    WAVEMANIPPAPI bool     isEmpty() const;
    WAVEMANIPPAPI Wave     peekWave() const;
    WAVEMANIPPAPI Wave     toWave() const;
};

} // namespace wm

#endif // !RECORDING_BUFFER_H
