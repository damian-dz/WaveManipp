#pragma once

/*!
 * \file WaveFileReader.hpp
 * \brief Bounded, random-access reader over a WAVE file's data chunk.
 */

#ifndef WAVE_FILE_READER_H
#define WAVE_FILE_READER_H

#include "utils.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace wm {
/*!
 * \brief Parses WAVE metadata and decodes bounded frame ranges without
 * loading the whole data chunk into memory.
 *
 * \details Unlike \ref Wave, which decodes and retains an entire file's
 * samples, WaveFileReader only parses the RIFF/fmt/data chunk headers at
 * open() time and then converts arbitrary frame ranges to normalized
 * interleaved floats on demand via readFrames(). Each instance owns its own
 * file handle and seek position -- independent readers over the same file
 * (e.g. one for playback, one for peak generation) never share a cursor or
 * lock. Move-only: copying would require either sharing a FILE* (unsafe,
 * shared seek position) or reopening the file (surprising), so neither is
 * offered.
 *
 * Supports the same formats as Wave: 8/16/24-bit PCM and 32-bit IEEE float,
 * RIFF (little-endian) or RIFX (big-endian), including WAVE_FORMAT_EXTENSIBLE
 * fmt chunks that resolve to one of those.
 */
class WAVEMANIPPAPI WaveFileReader
{
    static_assert(std::numeric_limits<float>::is_iec559, "float must be IEEE 754");
    static_assert(sizeof(float) == sizeof(uint32_t), "float must be 32 bits wide");

public:
    /*! \brief Creates a reader with no file open. */
    WaveFileReader() = default;
    /*! \brief Parses \p filename's header. \param filename Path to a RIFF or RIFX WAVE file. */
    explicit WaveFileReader(const char* filename);
    /*! \brief Parses \p filename's header. \param filename Path to a RIFF or RIFX WAVE file. */
    explicit WaveFileReader(const std::string& filename);

    WaveFileReader(const WaveFileReader&) = delete;
    WaveFileReader& operator=(const WaveFileReader&) = delete;
    /*! \brief Takes ownership of \p other's open file, leaving \p other closed. */
    WaveFileReader(WaveFileReader&& other) noexcept;
    /*! \brief Closes this reader's own file (if any), then takes ownership of \p other's. */
    WaveFileReader& operator=(WaveFileReader&& other) noexcept;
    /*! \brief Closes the underlying file, if open. */
    ~WaveFileReader();

    /*!
     * \brief Closes any currently open file, then parses \p filename's header.
     * \param filename Path to a RIFF or RIFX WAVE file.
     */
    void open(const char* filename);
    /*!
     * \brief Closes any currently open file, then parses \p filename's header.
     * \param filename Path to a RIFF or RIFX WAVE file.
     */
    void open(const std::string& filename);
    /*! \brief Closes the underlying file, if open. Safe to call when already closed. */
    void close();
    /*! \brief Returns true when a file is currently open. */
    bool isOpen() const noexcept;

    /*! \brief Returns the WAVE audio format code (1 = PCM, 3 = IEEE float). */
    uint16_t getAudioFormat() const noexcept;
    /*! \brief Returns the number of interleaved channels. */
    uint16_t getNumChannels() const noexcept;
    /*! \brief Returns the sample rate in Hz. */
    uint32_t getSampleRate() const noexcept;
    /*! \brief Returns the bit depth as stored in the file. */
    uint16_t getSampleBitDepth() const noexcept;
    /*! \brief Returns the number of audio frames in the data chunk. */
    uint32_t getNumFrames() const noexcept;
    /*! \brief Returns true when the file is RIFF (little-endian) rather than RIFX. */
    bool isLittleEndian() const noexcept;

    /*!
     * \brief Reads a bounded frame range, converted to normalized interleaved floats.
     * \param startFrame              Zero-based starting frame index.
     * \param frameCount              Maximum number of frames to read.
     * \param interleavedDestination  Destination buffer; must have room for at least
     *                                `frameCount * getNumChannels()` floats.
     * \details A `startFrame` at or past getNumFrames() reads zero frames; a read that
     * would run past the end (including one truncated by the file shrinking on disk
     * after open()) is truncated, not zero-filled -- matching AudioReader::readFrames's
     * contract (docs/AUDIO_STREAMING_PLAN.md).
     * \return the number of frames actually read (0 to \p frameCount).
     */
    uint32_t readFrames(uint64_t startFrame, uint32_t frameCount, float* interleavedDestination);

private:
    std::FILE* m_file = nullptr;
    long m_dataPosition = 0;
    bool m_isLittleEndian = true;
    uint16_t m_audioFormat = 0;
    uint16_t m_numChannels = 0;
    uint32_t m_sampleRate = 0;
    uint16_t m_bitsPerSample = 0;
    uint32_t m_numFrames = 0;
    std::vector<uint8_t> m_scratch;

    void moveFrom(WaveFileReader& other) noexcept;
};
}

#endif // !WAVE_FILE_READER_H
