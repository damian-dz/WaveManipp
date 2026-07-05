#include "WaveFileFormat.hpp"

#include "utils.hpp"

#include <cstring>
#include <vector>

namespace wm::detail {

namespace {

struct Descriptor {
    char id[4];
    uint32_t size;
};

struct FmtChunkBody {
    uint16_t audioFormat;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
};

} // namespace

bool isCpuBigEndian()
{
    union {
        uint16_t word;
        uint8_t bytes[2];
    } test = { 0x0A0B };
    return test.bytes[0] == 0x0A;
}

void changeBufferEndianness(uint8_t* buffer, uint32_t sampleSize, uint32_t bufferSize)
{
    for (uint32_t i = 0; i < bufferSize; i += sampleSize) {
        uint32_t j = i + sampleSize - 1;
        uint32_t k = i;
        while (j > k) {
            std::swap(buffer[j], buffer[k]);
            --j;
            ++k;
        }
    }
}

WaveFileInfo parseWaveFile(std::FILE* file)
{
    struct RiffHeader {
        Descriptor descriptor;
        char type[4];
    } riff{};

    if (std::fread(&riff, 1, sizeof(riff), file) != sizeof(riff)) {
        throwError("Wrong format or file corrupt.", "wm::detail::parseWaveFile");
    }

    bool isLittleEndian;
    if (std::memcmp(riff.descriptor.id, "RIFF", 4) == 0) {
        isLittleEndian = true;
    } else if (std::memcmp(riff.descriptor.id, "RIFX", 4) == 0) {
        isLittleEndian = false;
    } else {
        throwError("Wrong format or file corrupt.", "wm::detail::parseWaveFile");
    }
    if (std::memcmp(riff.type, "WAVE", 4) != 0) {
        throwError("Wrong format or file corrupt.", "wm::detail::parseWaveFile");
    }

    const bool isCpuLittleEndian = !isCpuBigEndian();
    const bool isEndiannessMismatched = isCpuLittleEndian != isLittleEndian;
    const auto decodeU32 = [isEndiannessMismatched](uint32_t value) {
        return isEndiannessMismatched ? reverseBytes<uint32_t>(value) : value;
    };
    const auto decodeU16 = [isEndiannessMismatched](uint16_t value) {
        return isEndiannessMismatched ? reverseBytes<uint16_t>(value) : value;
    };

    WaveFileInfo info;
    info.isLittleEndian = isLittleEndian;
    info.riffChunkSize = decodeU32(riff.descriptor.size);

    bool foundFmtChunk = false;
    bool foundDataChunk = false;
    long dataPosition = 0;
    uint32_t dataChunkSize = 0;
    FmtChunkBody fmt{};
    uint32_t fmtChunkSize = 0;
    const uint16_t waveFormatExtensible = 0xFFFE;

    while (!foundDataChunk || !foundFmtChunk) {
        Descriptor descriptor{};
        size_t bytesRead = std::fread(&descriptor, 1, sizeof(Descriptor), file);
        if (bytesRead == 0 && std::feof(file)) {
            break;
        }
        if (bytesRead != sizeof(Descriptor)) {
            throwError("Unexpected end of file while reading WAV chunks.", "wm::detail::parseWaveFile");
        }

        uint32_t chunkSize = decodeU32(descriptor.size);
        long chunkDataPosition = std::ftell(file);
        if (chunkDataPosition < 0) {
            throwError("Failed to inspect WAV chunk position.", "wm::detail::parseWaveFile");
        }

        if (std::memcmp(descriptor.id, "fmt ", 4) == 0) {
            if (chunkSize < 16) {
                throwError("Invalid WAV fmt chunk.", "wm::detail::parseWaveFile");
            }
            if (std::fread(&fmt, 1, sizeof(FmtChunkBody), file) != sizeof(FmtChunkBody)) {
                throwError("Unexpected end of file while reading WAV fmt chunk.", "wm::detail::parseWaveFile");
            }
            fmtChunkSize = chunkSize;
            fmt.audioFormat = decodeU16(fmt.audioFormat);
            fmt.numChannels = decodeU16(fmt.numChannels);
            fmt.sampleRate = decodeU32(fmt.sampleRate);
            fmt.byteRate = decodeU32(fmt.byteRate);
            fmt.blockAlign = decodeU16(fmt.blockAlign);
            fmt.bitsPerSample = decodeU16(fmt.bitsPerSample);

            std::vector<uint8_t> fmtExtension;
            uint32_t fmtExtensionSize = chunkSize - 16;
            if (fmtExtensionSize > 0) {
                fmtExtension.resize(fmtExtensionSize);
                if (std::fread(fmtExtension.data(), 1, fmtExtensionSize, file) != fmtExtensionSize) {
                    throwError("Unexpected end of file while reading WAV fmt extension.",
                               "wm::detail::parseWaveFile");
                }
            }

            if (fmt.audioFormat == waveFormatExtensible) {
                if (fmtExtension.size() < 24) {
                    throwError("Invalid WAVE_FORMAT_EXTENSIBLE fmt chunk.", "wm::detail::parseWaveFile");
                }
                const uint8_t* subFormat = fmtExtension.data() + 8;
                const bool isPcmSubFormat =
                    subFormat[0] == 0x01 && subFormat[1] == 0x00 && subFormat[2] == 0x00 && subFormat[3] == 0x00;
                const bool isFloatSubFormat =
                    subFormat[0] == 0x03 && subFormat[1] == 0x00 && subFormat[2] == 0x00 && subFormat[3] == 0x00;
                const bool hasWaveFormatGuidTail =
                    subFormat[4] == 0x00 && subFormat[5] == 0x00 && subFormat[6] == 0x10 && subFormat[7] == 0x00 &&
                    subFormat[8] == 0x80 && subFormat[9] == 0x00 && subFormat[10] == 0x00 && subFormat[11] == 0xAA &&
                    subFormat[12] == 0x00 && subFormat[13] == 0x38 && subFormat[14] == 0x9B && subFormat[15] == 0x71;
                if (!hasWaveFormatGuidTail || (!isPcmSubFormat && !isFloatSubFormat)) {
                    throwError("WAVE_FORMAT_EXTENSIBLE sub-format is not supported.", "wm::detail::parseWaveFile");
                }
                fmt.audioFormat = isPcmSubFormat ? 1 : 3;
            }

            foundFmtChunk = true;
        } else if (std::memcmp(descriptor.id, "data", 4) == 0) {
            dataChunkSize = chunkSize;
            dataPosition = std::ftell(file);
            foundDataChunk = true;
        }

        uint32_t paddingByte = chunkSize % 2;
        if (std::fseek(file, chunkDataPosition + chunkSize + paddingByte, SEEK_SET) != 0) {
            throwError("WAV chunk exceeds the file size.", "wm::detail::parseWaveFile");
        }
    }

    if (!foundFmtChunk || !foundDataChunk) {
        throwError("Required WAV fmt or data chunk is missing.", "wm::detail::parseWaveFile");
    }

    uint16_t bytesPerSample = fmt.bitsPerSample / 8;
    uint16_t frameSize = fmt.numChannels * bytesPerSample;

    if (fmt.numChannels == 0 || fmt.bitsPerSample == 0 || fmt.bitsPerSample % 8 != 0 || frameSize == 0) {
        throwError("Invalid WAV format parameters.", "wm::detail::parseWaveFile");
    }
    if ((fmt.audioFormat == 1 && fmt.bitsPerSample != 8 && fmt.bitsPerSample != 16 && fmt.bitsPerSample != 24) ||
        (fmt.audioFormat == 3 && fmt.bitsPerSample != 32) ||
        (fmt.audioFormat != 1 && fmt.audioFormat != 3)) {
        throwError("WAV audio format is not supported.", "wm::detail::parseWaveFile");
    }
    if (dataChunkSize % bytesPerSample != 0 || dataChunkSize % frameSize != 0) {
        throwError("WAV data chunk size does not align with the sample format.", "wm::detail::parseWaveFile");
    }

    info.fmtChunkSize = fmtChunkSize;
    info.audioFormat = fmt.audioFormat;
    info.numChannels = fmt.numChannels;
    info.sampleRate = fmt.sampleRate;
    info.byteRate = fmt.byteRate;
    info.blockAlign = fmt.blockAlign;
    info.bitsPerSample = fmt.bitsPerSample;
    info.dataChunkSize = dataChunkSize;
    info.dataOffset = static_cast<uint64_t>(dataPosition);
    return info;
}

} // namespace wm::detail
