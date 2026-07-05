#pragma once

/*!
 * \file WaveFileWriter.hpp
 * \brief Incremental, bounded-memory RIFF/RIFX WAVE writer.
 */

#ifndef WAVE_FILE_WRITER_H
#define WAVE_FILE_WRITER_H

#include "utils.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace wm {

/*! \brief WAVE sample encoding written by WaveFileWriter. */
enum class WaveFileEncoding : uint16_t {
    PcmInteger = 1,
    IeeeFloat = 3
};

/*! \brief Format of an incremental WAVE output stream. */
struct WaveFileWriterFormat {
    WaveFileEncoding encoding = WaveFileEncoding::PcmInteger;
    uint16_t numChannels = 2;
    uint32_t sampleRate = default_sample_rate;
    uint16_t bitsPerSample = 16;
    bool isLittleEndian = true;
};

/*!
 * \brief Writes normalized interleaved samples incrementally without retaining
 * previously written frames.
 *
 * \details Each successful writeFrames() call checkpoints the RIFF and data
 * sizes and flushes the stream, leaving a recoverable valid prefix if writing is
 * interrupted. finalize() appends the required odd-byte padding, writes the final
 * sizes, flushes, and closes the stream. Classic RIFF/RIFX size overflow is
 * rejected; RF64 is not yet supported.
 */
class WAVEMANIPPAPI WaveFileWriter
{
public:
    WaveFileWriter() = default;
    WaveFileWriter(const char* filename, const WaveFileWriterFormat& format);
    WaveFileWriter(const std::string& filename, const WaveFileWriterFormat& format);

    WaveFileWriter(const WaveFileWriter&) = delete;
    WaveFileWriter& operator=(const WaveFileWriter&) = delete;
    WaveFileWriter(WaveFileWriter&& other) noexcept;
    WaveFileWriter& operator=(WaveFileWriter&& other) noexcept;
    ~WaveFileWriter();

    void open(const char* filename, const WaveFileWriterFormat& format);
    void open(const std::string& filename, const WaveFileWriterFormat& format);

    /*! \brief Encodes, appends, checkpoints, and flushes interleaved float frames. */
    uint32_t writeFrames(const float* interleavedSource, uint32_t frameCount);
    /*! \brief Encodes, appends, checkpoints, and flushes interleaved double frames. */
    uint32_t writeFrames(const double* interleavedSource, uint32_t frameCount);

    /*! \brief Rewrites and flushes the current recoverable RIFF/data sizes. */
    void checkpoint();
    /*! \brief Writes final padding and sizes, flushes, and closes. Idempotent. */
    void finalize();
    /*! \brief Closes without adding final padding; the last successful checkpoint remains recoverable. */
    void close() noexcept;

    bool isOpen() const noexcept;
    const WaveFileWriterFormat& getFormat() const noexcept;
    uint64_t getNumFramesWritten() const noexcept;
    uint64_t getDataBytesWritten() const noexcept;

private:
    std::FILE* m_file = nullptr;
    WaveFileWriterFormat m_format{};
    uint32_t m_frameSize = 0;
    uint64_t m_numFrames = 0;
    uint64_t m_dataBytes = 0;
    bool m_failed = false;
    std::vector<uint8_t> m_scratch;

    template <typename Sample>
    uint32_t writeFramesImpl(const Sample* interleavedSource, uint32_t frameCount);

    void validateFormat(const WaveFileWriterFormat& format) const;
    void writeHeader(bool includeFinalPadding);
    void moveFrom(WaveFileWriter& other) noexcept;
};

} // namespace wm

#endif // !WAVE_FILE_WRITER_H
