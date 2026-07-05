#include "WaveFileReader.hpp"
#include "WaveFileFormat.hpp"

#include <algorithm>
#include <cstring>

namespace wm {

WaveFileReader::WaveFileReader(const char* filename)
{
    open(filename);
}

WaveFileReader::WaveFileReader(const std::string& filename)
{
    open(filename);
}

WaveFileReader::WaveFileReader(WaveFileReader&& other) noexcept
{
    moveFrom(other);
}

WaveFileReader& WaveFileReader::operator=(WaveFileReader&& other) noexcept
{
    if (this != &other) {
        close();
        moveFrom(other);
    }
    return *this;
}

void WaveFileReader::moveFrom(WaveFileReader& other) noexcept
{
    m_file = other.m_file;
    m_dataPosition = other.m_dataPosition;
    m_isLittleEndian = other.m_isLittleEndian;
    m_audioFormat = other.m_audioFormat;
    m_numChannels = other.m_numChannels;
    m_sampleRate = other.m_sampleRate;
    m_bitsPerSample = other.m_bitsPerSample;
    m_numFrames = other.m_numFrames;
    m_scratch = std::move(other.m_scratch);

    other.m_file = nullptr;
    other.m_dataPosition = 0;
    other.m_numFrames = 0;
}

WaveFileReader::~WaveFileReader()
{
    close();
}

void WaveFileReader::open(const std::string& filename)
{
    open(filename.c_str());
}

void WaveFileReader::open(const char* filename)
{
    close();

    detail::FileHandle file(std::fopen(filename, "rb"));
    if (!file) {
        throwError("Failed to open the file.", "WaveFileReader::open(const char*)");
    }

    // parseWaveFile() defines the supported-format/validation policy this class relies
    // on -- see Wave::open(), the other caller sharing it.
    const detail::WaveFileInfo info = detail::parseWaveFile(file.get());

    const uint16_t bytesPerSample = info.bitsPerSample / 8;
    const uint32_t frameSize = uint32_t(info.numChannels) * bytesPerSample;

    m_dataPosition = static_cast<long>(info.dataOffset);
    m_isLittleEndian = info.isLittleEndian;
    m_audioFormat = info.audioFormat;
    m_numChannels = info.numChannels;
    m_sampleRate = info.sampleRate;
    m_bitsPerSample = info.bitsPerSample;
    m_numFrames = info.dataChunkSize / frameSize;
    m_file = file.release();
}

void WaveFileReader::close()
{
    if (m_file) {
        std::fclose(m_file);
        m_file = nullptr;
    }
    m_dataPosition = 0;
    m_numFrames = 0;
}

bool WaveFileReader::isOpen() const noexcept
{
    return m_file != nullptr;
}

uint16_t WaveFileReader::getAudioFormat() const noexcept
{
    return m_audioFormat;
}

uint16_t WaveFileReader::getNumChannels() const noexcept
{
    return m_numChannels;
}

uint32_t WaveFileReader::getSampleRate() const noexcept
{
    return m_sampleRate;
}

uint16_t WaveFileReader::getSampleBitDepth() const noexcept
{
    return m_bitsPerSample;
}

uint32_t WaveFileReader::getNumFrames() const noexcept
{
    return m_numFrames;
}

bool WaveFileReader::isLittleEndian() const noexcept
{
    return m_isLittleEndian;
}

uint32_t WaveFileReader::readFrames(uint64_t startFrame, uint32_t frameCount, float* interleavedDestination)
{
    if (!m_file || frameCount == 0 || startFrame >= m_numFrames) {
        return 0;
    }

    const uint64_t available = uint64_t(m_numFrames) - startFrame;
    const uint32_t framesToRead = uint32_t(std::min<uint64_t>(available, frameCount));
    const uint16_t bytesPerSample = m_bitsPerSample / 8;
    const uint32_t frameSize = uint32_t(m_numChannels) * bytesPerSample;
    const uint64_t byteOffset = uint64_t(m_dataPosition) + startFrame * frameSize;

    if (std::fseek(m_file, static_cast<long>(byteOffset), SEEK_SET) != 0) {
        return 0;
    }

    const size_t bytesWanted = size_t(framesToRead) * frameSize;
    if (m_scratch.size() < bytesWanted) {
        m_scratch.resize(bytesWanted);
    }
    const size_t bytesRead = std::fread(m_scratch.data(), 1, bytesWanted, m_file);
    const uint32_t framesActuallyRead = uint32_t(bytesRead / frameSize);
    if (framesActuallyRead == 0) {
        return 0;
    }

    const bool isEndiannessMismatched = detail::isCpuBigEndian() == m_isLittleEndian;
    const uint32_t samplesToDecode = framesActuallyRead * uint32_t(m_numChannels);
    uint8_t* raw = m_scratch.data();

    if (m_bitsPerSample == 8) {
        for (uint32_t i = 0; i < samplesToDecode; ++i) {
            interleavedDestination[i] = raw[i] / 127.5f - 1.f;
        }
    } else if (m_bitsPerSample == 16) {
        if (isEndiannessMismatched) {
            detail::changeBufferEndianness(raw, 2, uint32_t(bytesRead));
        }
        const int16_t* src = reinterpret_cast<const int16_t*>(raw);
        for (uint32_t i = 0; i < samplesToDecode; ++i) {
            interleavedDestination[i] = src[i] / 32768.f;
        }
    } else if (m_bitsPerSample == 24) {
        if (isEndiannessMismatched) {
            detail::changeBufferEndianness(raw, 3, uint32_t(bytesRead));
        }
        for (uint32_t i = 0; i < samplesToDecode; ++i) {
            uint32_t j = i * 3;
            interleavedDestination[i] = (raw[j] << 8 | raw[j + 1] << 16 | raw[j + 2] << 24) / 2147483648.f;
        }
    } else if (m_bitsPerSample == 32) {
        if (isEndiannessMismatched) {
            detail::changeBufferEndianness(raw, 4, uint32_t(bytesRead));
        }
        std::memcpy(interleavedDestination, raw, size_t(samplesToDecode) * sizeof(float));
    }

    return framesActuallyRead;
}

} // namespace wm
