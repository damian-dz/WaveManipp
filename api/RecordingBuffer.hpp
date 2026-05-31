#pragma once

#ifndef RECORDING_BUFFER_H
#define RECORDING_BUFFER_H

#include "Wave.hpp"
#include <atomic>
#include <cstdint>
#include <vector>

namespace wm {

/*!
 * \brief Collects live recording samples and converts them to Wave objects.
 *
 * \details append() is intended for the audio thread only. peekWave() can be used
 * from another thread while recording to take a snapshot of frames published so far.
 * toWave() should be called after the recording stream has stopped.
 */
class RecordingBuffer {
    std::vector<float>    m_samples;
    std::atomic<uint32_t> m_frameCount { 0 };
    uint16_t              m_numChannels = 2;
    uint32_t              m_sampleRate  = 44100;
    uint16_t              m_bitDepth    = 16;

public:
    /*! \brief Creates an empty stereo, 44100 Hz, 16-bit recording buffer. */
    RecordingBuffer() = default;

    /*!
     * \brief Clears the buffer and sets the format for subsequent appended samples.
     * \param reserveFrames Optional frame capacity to reserve.
     */
    WAVEMANIPPAPI void     reset(uint16_t numChannels, uint32_t sampleRate, uint16_t bitDepth,
                                 uint32_t reserveFrames = 0);
    /*! \brief Appends interleaved normalized float samples for a number of frames. */
    WAVEMANIPPAPI void     append(const float* data, uint32_t numFrames);
    /*! \brief Returns the number of frames published by append(). */
    WAVEMANIPPAPI uint32_t numFrames() const;
    /*! \brief Returns true when no frames have been appended. */
    WAVEMANIPPAPI bool     isEmpty() const;
    /*! \brief Copies the currently published frames into a Wave snapshot. */
    WAVEMANIPPAPI Wave     peekWave() const;
    /*! \brief Converts the complete stopped recording buffer into a Wave. */
    WAVEMANIPPAPI Wave     toWave() const;
};

} // namespace wm

#endif // !RECORDING_BUFFER_H
