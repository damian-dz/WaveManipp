#include "WaveFileWriter.hpp"

#include "PcmCodec.hpp"
#include "WaveFileFormat.hpp"

#include <array>
#include <cstring>
#include <limits>
#include <utility>

namespace wm {

namespace {

constexpr uint64_t waveHeaderSize = 44;
constexpr uint64_t riffSizeWithoutData = 36;

void storeU16(uint8_t* destination, uint16_t value, bool littleEndian)
{
    destination[0] = uint8_t(littleEndian ? value : value >> 8);
    destination[1] = uint8_t(littleEndian ? value >> 8 : value);
}

void storeU32(uint8_t* destination, uint32_t value, bool littleEndian)
{
    destination[0] = uint8_t(littleEndian ? value : value >> 24);
    destination[1] = uint8_t(littleEndian ? value >> 8 : value >> 16);
    destination[2] = uint8_t(littleEndian ? value >> 16 : value >> 8);
    destination[3] = uint8_t(littleEndian ? value >> 24 : value);
}

} // namespace

WaveFileWriter::WaveFileWriter(const char* filename, const WaveFileWriterFormat& format)
{
    open(filename, format);
}

WaveFileWriter::WaveFileWriter(const std::string& filename, const WaveFileWriterFormat& format)
{
    open(filename, format);
}

WaveFileWriter::WaveFileWriter(WaveFileWriter&& other) noexcept
{
    moveFrom(other);
}

WaveFileWriter& WaveFileWriter::operator=(WaveFileWriter&& other) noexcept
{
    if (this != &other) {
        close();
        moveFrom(other);
    }
    return *this;
}

WaveFileWriter::~WaveFileWriter()
{
    close();
}

void WaveFileWriter::moveFrom(WaveFileWriter& other) noexcept
{
    m_file = other.m_file;
    m_format = other.m_format;
    m_frameSize = other.m_frameSize;
    m_numFrames = other.m_numFrames;
    m_dataBytes = other.m_dataBytes;
    m_failed = other.m_failed;
    m_scratch = std::move(other.m_scratch);

    other.m_file = nullptr;
    other.m_frameSize = 0;
    other.m_numFrames = 0;
    other.m_dataBytes = 0;
    other.m_failed = false;
}

void WaveFileWriter::validateFormat(const WaveFileWriterFormat& format) const
{
    const bool isPcm = format.encoding == WaveFileEncoding::PcmInteger;
    const bool isFloat = format.encoding == WaveFileEncoding::IeeeFloat;
    if ((!isPcm && !isFloat) || format.numChannels == 0 || format.sampleRate == 0) {
        throwError("Invalid WAVE writer format.", "WaveFileWriter::open");
    }
    if ((isPcm && format.bitsPerSample != 8 && format.bitsPerSample != 16 && format.bitsPerSample != 24) ||
        (isFloat && format.bitsPerSample != 32)) {
        throwError("WAVE writer encoding and bit depth are not supported.", "WaveFileWriter::open");
    }

    const uint64_t frameSize = uint64_t(format.numChannels) * format.bitsPerSample / 8;
    const uint64_t byteRate = uint64_t(format.sampleRate) * frameSize;
    if (frameSize == 0 || frameSize > std::numeric_limits<uint16_t>::max() ||
        byteRate > std::numeric_limits<uint32_t>::max()) {
        throwError("WAVE writer format exceeds RIFF field limits.", "WaveFileWriter::open");
    }
}

void WaveFileWriter::open(const char* filename, const WaveFileWriterFormat& format)
{
    close();
    validateFormat(format);
    if (!filename || filename[0] == '\0') {
        throwError("The WAVE output filename is empty.", "WaveFileWriter::open");
    }

    m_format = format;
    m_frameSize = uint32_t(format.numChannels) * format.bitsPerSample / 8;
    m_numFrames = 0;
    m_dataBytes = 0;
    m_failed = false;
    m_scratch.clear();
    m_file = std::fopen(filename, "wb+");
    if (!m_file) {
        throwError("Failed to create the file.", "WaveFileWriter::open");
    }

    try {
        writeHeader(false);
        if (!detail::seekFile(m_file, waveHeaderSize)) {
            m_failed = true;
            throwError("Failed to position the WAVE output stream.", "WaveFileWriter::open");
        }
    } catch (...) {
        close();
        throw;
    }
}

void WaveFileWriter::open(const std::string& filename, const WaveFileWriterFormat& format)
{
    open(filename.c_str(), format);
}

void WaveFileWriter::writeHeader(bool includeFinalPadding)
{
    const uint64_t padding = includeFinalPadding ? (m_dataBytes & 1u) : 0;
    const uint64_t riffSize = riffSizeWithoutData + m_dataBytes + padding;
    if (m_dataBytes > std::numeric_limits<uint32_t>::max() ||
        riffSize > std::numeric_limits<uint32_t>::max()) {
        m_failed = true;
        throwError("WAVE output exceeds the classic RIFF size limit.", "WaveFileWriter::writeHeader");
    }

    std::array<uint8_t, waveHeaderSize> header{};
    std::memcpy(header.data(), m_format.isLittleEndian ? "RIFF" : "RIFX", 4);
    storeU32(header.data() + 4, static_cast<uint32_t>(riffSize), m_format.isLittleEndian);
    std::memcpy(header.data() + 8, "WAVEfmt ", 8);
    storeU32(header.data() + 16, 16, m_format.isLittleEndian);
    storeU16(header.data() + 20, static_cast<uint16_t>(m_format.encoding), m_format.isLittleEndian);
    storeU16(header.data() + 22, m_format.numChannels, m_format.isLittleEndian);
    storeU32(header.data() + 24, m_format.sampleRate, m_format.isLittleEndian);
    storeU32(header.data() + 28, m_format.sampleRate * m_frameSize, m_format.isLittleEndian);
    storeU16(header.data() + 32, static_cast<uint16_t>(m_frameSize), m_format.isLittleEndian);
    storeU16(header.data() + 34, m_format.bitsPerSample, m_format.isLittleEndian);
    std::memcpy(header.data() + 36, "data", 4);
    storeU32(header.data() + 40, static_cast<uint32_t>(m_dataBytes), m_format.isLittleEndian);

    uint64_t previousPosition = 0;
    if (!detail::tellFile(m_file, previousPosition) || !detail::seekFile(m_file, 0) ||
        std::fwrite(header.data(), 1, header.size(), m_file) != header.size() || std::fflush(m_file) != 0 ||
        !detail::seekFile(m_file, previousPosition)) {
        m_failed = true;
        throwError("Failed to checkpoint the WAVE header.", "WaveFileWriter::writeHeader");
    }
}

template <typename Sample>
uint32_t WaveFileWriter::writeFramesImpl(const Sample* interleavedSource, uint32_t frameCount)
{
    if (!m_file) {
        throwError("No WAVE output file is open.", "WaveFileWriter::writeFrames");
    }
    if (m_failed) {
        throwError("The WAVE output stream is in a failed state.", "WaveFileWriter::writeFrames");
    }
    if (frameCount == 0) {
        return 0;
    }
    if (!interleavedSource) {
        throwError("The interleaved source buffer is null.", "WaveFileWriter::writeFrames");
    }

    const uint64_t bytesToWrite = uint64_t(frameCount) * m_frameSize;
    const uint64_t newDataBytes = m_dataBytes + bytesToWrite;
    const uint64_t finalRiffSize = riffSizeWithoutData + newDataBytes + (newDataBytes & 1u);
    if (newDataBytes < m_dataBytes || finalRiffSize > std::numeric_limits<uint32_t>::max() ||
        bytesToWrite > std::numeric_limits<size_t>::max()) {
        throwError("WAVE output exceeds the classic RIFF size limit.", "WaveFileWriter::writeFrames");
    }

    m_scratch.resize(static_cast<size_t>(bytesToWrite));
    const size_t sampleCount = size_t(frameCount) * m_format.numChannels;
    if (m_format.bitsPerSample == 8) {
        for (size_t i = 0; i < sampleCount; ++i) {
            m_scratch[i] = detail::encodePcm8Sample(interleavedSource[i]);
        }
    } else if (m_format.bitsPerSample == 16) {
        for (size_t i = 0; i < sampleCount; ++i) {
            const int32_t sample = detail::quantizeSignedPcm<Sample, 16>(interleavedSource[i]);
            storeU16(m_scratch.data() + i * 2, static_cast<uint16_t>(sample), m_format.isLittleEndian);
        }
    } else if (m_format.bitsPerSample == 24) {
        for (size_t i = 0; i < sampleCount; ++i) {
            detail::encodePcm24Sample(interleavedSource[i], m_format.isLittleEndian, m_scratch.data() + i * 3);
        }
    } else {
        for (size_t i = 0; i < sampleCount; ++i) {
            const float sample = static_cast<float>(interleavedSource[i]);
            uint32_t bits = 0;
            std::memcpy(&bits, &sample, sizeof(bits));
            storeU32(m_scratch.data() + i * 4, bits, m_format.isLittleEndian);
        }
    }

    if (std::fwrite(m_scratch.data(), 1, m_scratch.size(), m_file) != m_scratch.size()) {
        std::clearerr(m_file);
        if ((m_dataBytes & 1u) != 0 && detail::seekFile(m_file, waveHeaderSize + m_dataBytes)) {
            std::fputc(0, m_file);
            std::fflush(m_file);
        }
        m_failed = true;
        throwError("Failed to write WAVE sample data.", "WaveFileWriter::writeFrames");
    }

    m_dataBytes = newDataBytes;
    m_numFrames += frameCount;
    checkpoint();
    return frameCount;
}

uint32_t WaveFileWriter::writeFrames(const float* interleavedSource, uint32_t frameCount)
{
    return writeFramesImpl(interleavedSource, frameCount);
}

uint32_t WaveFileWriter::writeFrames(const double* interleavedSource, uint32_t frameCount)
{
    return writeFramesImpl(interleavedSource, frameCount);
}

void WaveFileWriter::checkpoint()
{
    if (!m_file) {
        throwError("No WAVE output file is open.", "WaveFileWriter::checkpoint");
    }
    if (m_failed) {
        throwError("The WAVE output stream is in a failed state.", "WaveFileWriter::checkpoint");
    }

    const uint64_t dataEnd = waveHeaderSize + m_dataBytes;
    if (!detail::seekFile(m_file, dataEnd)) {
        m_failed = true;
        throwError("Failed to position the WAVE output stream.", "WaveFileWriter::checkpoint");
    }
    if ((m_dataBytes & 1u) != 0 && std::fputc(0, m_file) == EOF) {
        m_failed = true;
        throwError("Failed to write WAVE padding.", "WaveFileWriter::checkpoint");
    }

    writeHeader(true);
    if (!detail::seekFile(m_file, dataEnd)) {
        m_failed = true;
        throwError("Failed to restore the WAVE data position.", "WaveFileWriter::checkpoint");
    }
}

void WaveFileWriter::finalize()
{
    if (!m_file) {
        return;
    }
    if (m_failed) {
        close();
        throwError("Cannot finalize a failed WAVE output stream.", "WaveFileWriter::finalize");
    }

    try {
        checkpoint();
    } catch (...) {
        close();
        throw;
    }

    std::FILE* file = m_file;
    m_file = nullptr;
    if (std::fclose(file) != 0) {
        m_failed = true;
        throwError("Failed to close the finalized WAVE output.", "WaveFileWriter::finalize");
    }
}

void WaveFileWriter::close() noexcept
{
    if (m_file) {
        std::fclose(m_file);
        m_file = nullptr;
    }
}

bool WaveFileWriter::isOpen() const noexcept
{
    return m_file != nullptr;
}

const WaveFileWriterFormat& WaveFileWriter::getFormat() const noexcept
{
    return m_format;
}

uint64_t WaveFileWriter::getNumFramesWritten() const noexcept
{
    return m_numFrames;
}

uint64_t WaveFileWriter::getDataBytesWritten() const noexcept
{
    return m_dataBytes;
}

} // namespace wm
