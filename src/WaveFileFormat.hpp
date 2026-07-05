#pragma once

/*!
 * \file WaveFileFormat.hpp
 * \brief Private RIFF/RIFX WAVE header parser shared by Wave and WaveFileReader.
 *
 * \details Not part of the public API (lives in src/, not api/): this is the one
 * place that defines which WAVE encodings are supported and how a header is
 * validated. Wave::open() and WaveFileReader::open() both delegate to
 * parseWaveFile() instead of each maintaining their own copy of that policy.
 */

#include <algorithm>
#include <cstdint>
#include <cstdio>

namespace wm::detail {

//! Decoded RIFF/fmt/data chunk metadata for one WAVE file. audioFormat is
//! already resolved from WAVE_FORMAT_EXTENSIBLE to a plain 1 (PCM) or 3
//! (IEEE float) if the fmt chunk used it.
struct WaveFileInfo {
    bool isLittleEndian = true;
    uint32_t riffChunkSize = 0;
    uint32_t fmtChunkSize = 0;
    uint16_t audioFormat = 0;
    uint16_t numChannels = 0;
    uint32_t sampleRate = 0;
    uint32_t byteRate = 0;
    uint16_t blockAlign = 0;
    uint16_t bitsPerSample = 0;
    uint32_t dataChunkSize = 0;
    uint64_t dataOffset = 0;
};

/*!
 * \brief Parses and validates a RIFF/RIFX WAVE file's header.
 * \param file An open, readable file positioned at the very start of the file.
 * \details Walks every chunk until both "fmt " and "data" are found, skipping
 * (and padding-byte-aligning past) anything else -- so an unrecognized chunk
 * anywhere in the layout is tolerated. Validates the same constraints
 * Wave::open() has always enforced: a supported bit depth for the resolved
 * audio format (8/16/24-bit PCM or 32-bit float), a data chunk size aligned
 * to the sample/frame size, and both required chunks present.
 *
 * On any malformed, truncated, or unsupported input this throws (via
 * wm::throwError) rather than closing \p file itself -- file ownership is
 * the caller's responsibility (see FileHandle), so a thrown exception during
 * parsing still needs the caller's RAII to close the file during unwinding.
 *
 * Leaves \p file at an unspecified position; callers that need to read
 * sample data should seek explicitly using the returned dataOffset.
 */
WaveFileInfo parseWaveFile(std::FILE* file);

/*!
 * \brief Seeks to an absolute byte offset without narrowing to 32-bit long.
 * \return true on success, false when the offset is not representable or the
 * underlying seek fails.
 */
bool seekFile(std::FILE* file, uint64_t offset) noexcept;

/*!
 * \brief Reports the current absolute byte offset without narrowing it.
 * \param file   Open file whose position is queried.
 * \param offset Receives the non-negative byte offset on success.
 * \return true on success.
 */
bool tellFile(std::FILE* file, uint64_t& offset) noexcept;

//! Returns true when the running CPU is big-endian.
bool isCpuBigEndian();

//! Reverses the byte order of an integral value.
template <typename T>
T reverseBytes(T val)
{
    char* pVal = reinterpret_cast<char*>(&val);
    T res = 0;
    char* pRes = reinterpret_cast<char*>(&res);
    std::reverse_copy(pVal, pVal + sizeof(T), pRes);
    return res;
}

//! Reverses the byte order of every `sampleSize`-byte sample within `buffer`.
void changeBufferEndianness(uint8_t* buffer, uint32_t sampleSize, uint32_t bufferSize);

/*!
 * \brief Minimal RAII owner for a std::FILE*.
 * \details Exists so parseWaveFile() (and any other fallible step between
 * fopen() and a class taking permanent ownership of the handle) can simply
 * throw on error without itself deciding whether to close the file --
 * closing happens here during stack unwinding, or via release() once the
 * caller is ready to keep the handle longer-term (e.g. WaveFileReader).
 */
class FileHandle {
public:
    FileHandle() = default;
    explicit FileHandle(std::FILE* file) : m_file(file) {}
    FileHandle(const FileHandle&) = delete;
    FileHandle& operator=(const FileHandle&) = delete;
    FileHandle(FileHandle&& other) noexcept : m_file(other.release()) {}
    FileHandle& operator=(FileHandle&& other) noexcept
    {
        if (this != &other) {
            reset();
            m_file = other.release();
        }
        return *this;
    }
    ~FileHandle() { reset(); }

    std::FILE* get() const noexcept { return m_file; }
    explicit operator bool() const noexcept { return m_file != nullptr; }

    //! Relinquishes ownership without closing the file.
    std::FILE* release() noexcept
    {
        std::FILE* file = m_file;
        m_file = nullptr;
        return file;
    }

    //! Closes the current file (if any) and takes ownership of `file` (if any).
    void reset(std::FILE* file = nullptr) noexcept
    {
        if (m_file) {
            std::fclose(m_file);
        }
        m_file = file;
    }

private:
    std::FILE* m_file = nullptr;
};

} // namespace wm::detail
