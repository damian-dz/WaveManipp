#include "Wave.hpp"
#include "Dsp.hpp"
#include <cmath>

namespace wm {

Wave::Wave() :
    m_isLittleEndian(true),
    m_numSamples(0)
{
    zeroInitHeader();
}

Wave::Wave(const char* filename) :
    m_isLittleEndian(true),
    m_numSamples(0)
{
    zeroInitHeader();
    open(filename);
}

Wave::Wave(const std::string& filename) :
    m_isLittleEndian(true),
    m_numSamples(0)
{
    zeroInitHeader();
    open(filename);
}

Wave::Wave(uint32_t numFrames, uint16_t numChannels, uint16_t bitDepth, uint32_t sampleRate) :
    m_isLittleEndian(true),
    m_numSamples(numChannels * numFrames)
{
    zeroInitHeader();
    setFourCharacterCodes();
    m_waveProperties.setNumChannels(numChannels);
    reserveMemory(m_numSamples);
    setSampleRate(sampleRate);
    setSampleBitDepth(bitDepth);
    m_waveProperties.setRiffChunkSize(m_waveProperties.getDataChunkSize() + 36);
    m_waveProperties.setFmtChunkSize(16);
}

Wave::Wave(const Wave& other) :
    m_header(other.m_header),
    m_dataSubChunk(other.m_dataSubChunk),
    m_isLittleEndian(other.m_isLittleEndian),
    m_numSamples(other.m_numSamples),
    m_buffer(other.m_buffer),
    m_waveProperties(other.m_waveProperties)
{}

Wave::~Wave() = default;

void Wave::setFourCharacterCodes()
{
    if (m_isLittleEndian) {
        std::memcpy(m_header.riff.descriptor.id, "RIFF", 4);
    } else {
        std::memcpy(m_header.riff.descriptor.id, "RIFX", 4);
    }
    std::memcpy(m_header.riff.type, "WAVE", 4);
    std::memcpy(m_header.wave.descriptor.id, "fmt ", 4);
    std::memcpy(m_dataSubChunk.descriptor.id, "data", 4);
}

void Wave::generateHeader()
{
    setFourCharacterCodes();

    bool isCpuLittleEndian = !isCpuBigEndian();
    bool isEndiannessMismatched = isCpuLittleEndian != m_isLittleEndian;

    m_header.riff.descriptor.size = isEndiannessMismatched ?
        reverseBytes<uint32_t>(m_waveProperties.getRiffChunkSize()) : m_waveProperties.getRiffChunkSize();
    m_header.wave.descriptor.size = isEndiannessMismatched ?
        reverseBytes<uint32_t>(m_waveProperties.getFmtChunkSize()) : m_waveProperties.getFmtChunkSize();
    m_header.wave.audioFormat = isEndiannessMismatched ?
        reverseBytes<uint16_t>(m_waveProperties.getAudioFormat()) : m_waveProperties.getAudioFormat();
    m_header.wave.numChannels = isEndiannessMismatched ?
        reverseBytes<uint16_t>(m_waveProperties.getNumChannels()) : m_waveProperties.getNumChannels();
    m_header.wave.sampleRate = isEndiannessMismatched ?
        reverseBytes<uint32_t>(m_waveProperties.getSamplingFrequency()) : m_waveProperties.getSamplingFrequency();
    m_header.wave.byteRate = isEndiannessMismatched ?
        reverseBytes<uint32_t>(m_waveProperties.getNumBytesPerSecond()) : m_waveProperties.getNumBytesPerSecond();
    m_header.wave.blockAlign = isEndiannessMismatched ?
        reverseBytes<uint16_t>(m_waveProperties.getBlockAlign()) : m_waveProperties.getBlockAlign();
    m_header.wave.bitsPerSample = isEndiannessMismatched ?
        reverseBytes<uint16_t>(m_waveProperties.getNumBitsPerSample()) : m_waveProperties.getNumBitsPerSample();

    m_dataSubChunk.descriptor.size = isEndiannessMismatched ?
        reverseBytes<uint32_t>(m_waveProperties.getDataChunkSize()) : m_waveProperties.getDataChunkSize();
}

void Wave::changeBufferEndianness(uint8_t* buffer, uint32_t sampleSize, uint32_t bufferSize)
{
    for (uint32_t i = 0; i < bufferSize; i += sampleSize) {
        uint32_t j = i + sampleSize - 1;
        uint32_t k = i;
        while (j > k) {
            uint8_t u8 = buffer[j];
            buffer[j] = buffer[k];
            buffer[k] = u8;
            --j;
            ++k;
        }
    }
}

bool Wave::isCpuBigEndian() const
{
    union {
        uint16_t word;
        uint8_t bytes[2];
    } TestStruct = { 0x0A0B };
    return TestStruct.bytes[0] == 0x0A;
}

void Wave::setWaveProperties(bool onlyDataChunk)
{
    bool isCpuLittleEndian = !isCpuBigEndian();
    bool isEndiannessMismatched = isCpuLittleEndian != m_isLittleEndian;
    if (!onlyDataChunk) {
        m_waveProperties.setRiffChunkSize(isEndiannessMismatched ?
            reverseBytes<uint32_t>(m_header.riff.descriptor.size) : m_header.riff.descriptor.size);
        m_waveProperties.setFmtChunkSize(isEndiannessMismatched ?
            reverseBytes<uint32_t>(m_header.wave.descriptor.size) : m_header.wave.descriptor.size);
        m_waveProperties.setAudioFormat(isEndiannessMismatched ?
            reverseBytes<uint16_t>(m_header.wave.audioFormat) : m_header.wave.audioFormat);
        m_waveProperties.setNumChannels(isEndiannessMismatched ?
            reverseBytes<uint16_t>(m_header.wave.numChannels) : m_header.wave.numChannels);
        m_waveProperties.setSamplingFrequency(isEndiannessMismatched ?
            reverseBytes<uint32_t>(m_header.wave.sampleRate) : m_header.wave.sampleRate);
        m_waveProperties.setNumBytesPerSecond(isEndiannessMismatched ?
            reverseBytes<uint32_t>(m_header.wave.byteRate) : m_header.wave.byteRate);
        m_waveProperties.setBlockAlign(isEndiannessMismatched ?
            reverseBytes<uint16_t>(m_header.wave.blockAlign) : m_header.wave.blockAlign);
        m_waveProperties.setNumBitsPerSample(isEndiannessMismatched ?
            reverseBytes<uint16_t>(m_header.wave.bitsPerSample) : m_header.wave.bitsPerSample);
    } else {
        m_waveProperties.setDataChunkSize(isEndiannessMismatched ?
            reverseBytes<uint32_t>(m_dataSubChunk.descriptor.size) : m_dataSubChunk.descriptor.size);
    }
}

void Wave::detach()
{
    if (m_buffer && m_buffer.use_count() > 1) {
        auto copy = std::shared_ptr<float[]>(new float[m_numSamples]);
        std::memcpy(copy.get(), m_buffer.get(), m_numSamples * sizeof(float));
        m_buffer = std::move(copy);
    }
}

void Wave::reserveMemory(uint32_t numSamples, bool zeroInit)
{
    if (numSamples == 0) {
        m_buffer.reset();
        return;
    }
    m_buffer = std::shared_ptr<float[]>(new float[numSamples]);
    if (zeroInit) {
        std::memset(m_buffer.get(), 0, numSamples * sizeof(float));
    }
}

void Wave::resizeMemory(uint32_t numSamples, bool zeroInit)
{
    if (numSamples == m_numSamples) {
        throwError("The number of samples remains the same. Resize unnecessary.",
                   "Wave::resizeMemory(uint32_t, bool)");
    }
    uint32_t count = std::min(m_numSamples, numSamples);
    m_numSamples = numSamples;
    m_waveProperties.setDataChunkSize(numSamples * m_waveProperties.getNumBitsPerSample() / 8);
    m_waveProperties.setRiffChunkSize(m_waveProperties.getDataChunkSize() + 36);
    if (numSamples == 0) {
        m_buffer.reset();
    } else {
        auto newBuffer = std::shared_ptr<float[]>(new float[numSamples]);
        if (m_buffer && count > 0) {
            std::memcpy(newBuffer.get(), m_buffer.get(), count * sizeof(float));
        }
        if (zeroInit && numSamples > count) {
            std::memset(newBuffer.get() + count, 0, (numSamples - count) * sizeof(float));
        }
        m_buffer = std::move(newBuffer);
    }
}

void Wave::copySamples(const float* source, float* destination, uint32_t count, uint32_t srcOffset,
                       uint32_t destOffset)
{
    std::memcpy(&destination[destOffset], &source[srcOffset], count * sizeof(float));
}

bool Wave::peekForId(const std::string& id, std::FILE* pFile)
{
    bool res = true;
    int pos = std::ftell(pFile);
    for (size_t i = 0; i < id.size(); ++i) {
        int c = std::fgetc(pFile);
        if (id[i] != char(c)) {
            res = false;
            break;
        }
    }
    std::fseek(pFile, pos, 0);
    return res;
}

void Wave::open(const char* filename, uint32_t bufferSize)
{
    std::FILE* pFile = std::fopen(filename, "rb");
    if (!pFile) {
        throwError("Failed to open the file.",
                   "Wave::open(const std::string&, uint32_t)");
    }
    if (std::fread(&m_header.riff, 1, sizeof(RIFFChunk), pFile) != sizeof(RIFFChunk)) {
        std::fclose(pFile);
        throwError("Wrong format or file corrupt.",
                   "Wave::open(const std::string&, uint32_t)");
    }

    if (std::memcmp(m_header.riff.descriptor.id, "RIFF", 4) == 0) {
        m_isLittleEndian = true;
    } else if (std::memcmp(m_header.riff.descriptor.id, "RIFX", 4) == 0) {
        m_isLittleEndian = false;
    } else {
        std::fclose(pFile);
        throwError("Wrong format or file corrupt.",
                   "Wave::open(const std::string&, uint32_t)");
    }

    if (std::memcmp(m_header.riff.type, "WAVE", 4) != 0) {
        std::fclose(pFile);
        throwError("Wrong format or file corrupt.",
                   "Wave::open(const std::string&, uint32_t)");
    }

    const bool isCpuLittleEndian = !isCpuBigEndian();
    const bool isEndiannessMismatched = isCpuLittleEndian != m_isLittleEndian;
    const auto decodeU32 = [isEndiannessMismatched, this](uint32_t value) {
        return isEndiannessMismatched ? reverseBytes<uint32_t>(value) : value;
    };

    m_waveProperties.setRiffChunkSize(decodeU32(m_header.riff.descriptor.size));

    bool foundFmtChunk = false;
    bool foundDataChunk = false;
    long dataPosition = 0;
    const uint16_t waveFormatExtensible = 0xFFFE;

    while (!foundDataChunk || !foundFmtChunk) {
        Descriptor descriptor {};
        size_t bytesRead = std::fread(&descriptor, 1, sizeof(Descriptor), pFile);
        if (bytesRead == 0 && std::feof(pFile)) {
            break;
        }
        if (bytesRead != sizeof(Descriptor)) {
            std::fclose(pFile);
            throwError("Unexpected end of file while reading WAV chunks.",
                       "Wave::open(const std::string&, uint32_t)");
        }

        uint32_t chunkSize = decodeU32(descriptor.size);
        long chunkDataPosition = std::ftell(pFile);
        if (chunkDataPosition < 0) {
            std::fclose(pFile);
            throwError("Failed to inspect WAV chunk position.",
                       "Wave::open(const std::string&, uint32_t)");
        }

        if (std::memcmp(descriptor.id, "fmt ", 4) == 0) {
            if (chunkSize < 16) {
                std::fclose(pFile);
                throwError("Invalid WAV fmt chunk.",
                           "Wave::open(const std::string&, uint32_t)");
            }

            m_header.wave.descriptor = descriptor;
            if (std::fread(&m_header.wave.audioFormat, 1, sizeof(FmtSubChunk) - sizeof(Descriptor), pFile) !=
                sizeof(FmtSubChunk) - sizeof(Descriptor)) {
                std::fclose(pFile);
                throwError("Unexpected end of file while reading WAV fmt chunk.",
                           "Wave::open(const std::string&, uint32_t)");
            }

            std::vector<uint8_t> fmtExtension;
            uint32_t fmtExtensionSize = chunkSize - 16;
            if (fmtExtensionSize > 0) {
                fmtExtension.resize(fmtExtensionSize);
                if (std::fread(fmtExtension.data(), 1, fmtExtensionSize, pFile) != fmtExtensionSize) {
                    std::fclose(pFile);
                    throwError("Unexpected end of file while reading WAV fmt extension.",
                               "Wave::open(const std::string&, uint32_t)");
                }
            }

            setWaveProperties(false);

            if (m_waveProperties.getAudioFormat() == waveFormatExtensible) {
                if (fmtExtension.size() < 24) {
                    std::fclose(pFile);
                    throwError("Invalid WAVE_FORMAT_EXTENSIBLE fmt chunk.",
                               "Wave::open(const std::string&, uint32_t)");
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
                    std::fclose(pFile);
                    throwError("WAVE_FORMAT_EXTENSIBLE sub-format is not supported.",
                               "Wave::open(const std::string&, uint32_t)");
                }

                m_waveProperties.setAudioFormat(isPcmSubFormat ? 1 : 3);
            }

            foundFmtChunk = true;
        } else if (std::memcmp(descriptor.id, "data", 4) == 0) {
            m_dataSubChunk.descriptor = descriptor;
            m_waveProperties.setDataChunkSize(chunkSize);
            dataPosition = std::ftell(pFile);
            foundDataChunk = true;
        }

        uint32_t paddingByte = chunkSize % 2;
        if (std::fseek(pFile, chunkDataPosition + chunkSize + paddingByte, SEEK_SET) != 0) {
            std::fclose(pFile);
            throwError("WAV chunk exceeds the file size.",
                       "Wave::open(const std::string&, uint32_t)");
        }
    }

    if (!foundFmtChunk || !foundDataChunk) {
        std::fclose(pFile);
        throwError("Required WAV fmt or data chunk is missing.",
                   "Wave::open(const std::string&, uint32_t)");
    }

    uint16_t audioFormat = m_waveProperties.getAudioFormat();
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint16_t sampleBitDepth = m_waveProperties.getNumBitsPerSample();
    uint16_t bytesPerSample = sampleBitDepth / 8;
    uint16_t frameSize = numChannels * bytesPerSample;
    uint32_t dataChunkSize = m_waveProperties.getDataChunkSize();

    if (numChannels == 0 || sampleBitDepth == 0 || sampleBitDepth % 8 != 0 || frameSize == 0) {
        std::fclose(pFile);
        throwError("Invalid WAV format parameters.",
                   "Wave::open(const std::string&, uint32_t)");
    }
    if ((audioFormat == 1 && sampleBitDepth != 8 && sampleBitDepth != 16 && sampleBitDepth != 24) ||
        (audioFormat == 3 && sampleBitDepth != 32) ||
        (audioFormat != 1 && audioFormat != 3)) {
        std::fclose(pFile);
        throwError("WAV audio format is not supported.",
                   "Wave::open(const std::string&, uint32_t)");
    }
    if (dataChunkSize % bytesPerSample != 0 || dataChunkSize % frameSize != 0) {
        std::fclose(pFile);
        throwError("WAV data chunk size does not align with the sample format.",
                   "Wave::open(const std::string&, uint32_t)");
    }

    m_numSamples = dataChunkSize / bytesPerSample;
    reserveMemory(m_numSamples);
    if (m_buffer) {
        if (std::fseek(pFile, dataPosition, SEEK_SET) != 0) {
            std::fclose(pFile);
            throwError("Failed to seek to WAV data chunk.",
                       "Wave::open(const std::string&, uint32_t)");
        }
        readData(pFile, bufferSize);
    }
    std::fclose(pFile);
}

void Wave::open(const std::string& filename, uint32_t bufferSize)
{
    open(filename.c_str(), bufferSize);
}

void Wave::readData(std::FILE* file, uint32_t bufferSize)
{
    uint8_t* pBuffer = reinterpret_cast<uint8_t*>(std::malloc(bufferSize));
    if (pBuffer == nullptr) {
        throwError("Failed to create a buffer.",
                   "Wave::readData(std::FILE*, uint32_t)");
    }

    bool isEndiannessMismatched = isCpuBigEndian() == m_isLittleEndian;

    uint16_t sampleBitDepth = m_waveProperties.getNumBitsPerSample();
    uint32_t dataChunkSize = m_waveProperties.getDataChunkSize();
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint16_t frameSize = numChannels * sampleBitDepth / 8;
    if (sampleBitDepth == 0 || sampleBitDepth % 8 != 0 || numChannels == 0 || frameSize == 0) {
        std::free(pBuffer);
        throwError("Invalid WAV format parameters.",
                   "Wave::readData(std::FILE*, uint32_t)");
    }
    uint32_t numSamples = dataChunkSize / (uint32_t(sampleBitDepth) / 8);
    if (bufferSize % frameSize != 0) {
        std::free(pBuffer);
        throwError("The buffer size must be a multiple of the frame size.",
                   "Wave::readData(std::FILE*, uint32_t)");
    }

    float* dst = m_buffer.get();
    uint32_t offset = 0;
    if (sampleBitDepth == 8) {
        while (offset < numSamples) {
            uint32_t numBytesToRead = std::min(bufferSize, dataChunkSize - offset);
            std::fread(pBuffer, 1, numBytesToRead, file);
            uint32_t numBufferFrames = numBytesToRead / frameSize;
            for (uint32_t y = 0; y < numBufferFrames; ++y) {
                for (uint16_t x = 0; x < numChannels; ++x) {
                    uint32_t i = y * numChannels + x;
                    dst[offset + i] = pBuffer[i] / 127.5f - 1.f;
                }
            }
            offset += bufferSize;
        }
    } else if (sampleBitDepth == 16) {
        while (offset < numSamples) {
            uint32_t numBytesToRead = std::min(bufferSize, dataChunkSize - 2 * offset);
            std::fread(pBuffer, 1, numBytesToRead, file);
            if (isEndiannessMismatched) {
                changeBufferEndianness(pBuffer, 2, numBytesToRead);
            }
            uint32_t numBufferFrames = numBytesToRead / frameSize;
            int16_t* pI16Buffer = reinterpret_cast<int16_t*>(pBuffer);
            for (uint32_t y = 0; y < numBufferFrames; ++y) {
                for (uint16_t x = 0; x < numChannels; ++x) {
                    uint32_t i = y * numChannels + x;
                    dst[offset + i] = pI16Buffer[i] / 32768.f;
                }
            }
            offset += bufferSize / 2;
        }
    } else if (sampleBitDepth == 24) {
        while (offset < numSamples) {
            uint32_t numBytesToRead = std::min(bufferSize, dataChunkSize - 3 * offset);
            std::fread(pBuffer, 1, numBytesToRead, file);
            if (isEndiannessMismatched) {
                changeBufferEndianness(pBuffer, 3, numBytesToRead);
            }
            uint32_t numBufferFrames = numBytesToRead / frameSize;
            for (uint32_t y = 0; y < numBufferFrames; ++y) {
                for (uint16_t x = 0; x < numChannels; ++x) {
                    uint32_t i = 3 * (y * numChannels + x);
                    dst[offset + i / 3] =
                        (pBuffer[i] << 8 | pBuffer[i + 1] << 16 | pBuffer[i + 2] << 24) / 2147483648.f;
                }
            }
            offset += bufferSize / 3;
        }
    } else if (sampleBitDepth == 32) {
        while (offset < numSamples) {
            uint32_t numBytesToRead = std::min(bufferSize, dataChunkSize - 4 * offset);
            std::fread(pBuffer, 1, numBytesToRead, file);
            if (isEndiannessMismatched) {
                changeBufferEndianness(pBuffer, 4, numBytesToRead);
            }
            uint32_t numBufferFrames = numBytesToRead / frameSize;
            float* pF32Buffer = reinterpret_cast<float*>(pBuffer);
            copySamples(pF32Buffer, dst, numBufferFrames * numChannels, 0, offset);
            offset += bufferSize / 4;
        }
    } else {
        std::free(pBuffer);
        throwError("Bit depth not supported.",
                   "Wave::readData(std::FILE*, size_t)");
    }
    std::free(pBuffer);
}

const float* Wave::constAudioData() const
{
    return m_buffer.get();
}

float* Wave::audioData()
{
    detach();
    return m_buffer.get();
}

std::vector<float> Wave::getBuffer(uint32_t offset, uint32_t sampleCount, int channel) const
{
    std::vector<float> buffer(sampleCount);
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint32_t absoluteOffset = numChannels * offset;
    for (uint32_t i = channel; i < sampleCount * numChannels; i += numChannels) {
        buffer[i / numChannels] = m_buffer[absoluteOffset + i];
    }
    return buffer;
}

std::vector<float> Wave::getSqueezedBuffer(uint32_t offset, uint32_t squeezedSampleCount, float squeezeFactor,
                                           bool absolute, int channel, bool multiThreaded) const
{
    std::vector<float> squeezedBuffer(squeezedSampleCount);
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint32_t absoluteOffset = numChannels * offset;
    #pragma omp parallel for if (multiThreaded)
    for (int32_t i = 0; i < int32_t(squeezedSampleCount); ++i) {
        int32_t current = int32_t(roundf(i * squeezeFactor)) * numChannels;
        int32_t next = int32_t(roundf((i + 1) * squeezeFactor)) * numChannels;
        float sum = 0.f;
        if (!absolute) {
            for (int32_t j = current; j < next; j += numChannels) {
                sum += m_buffer[absoluteOffset + j + channel];
            }
        } else {
            for (int32_t j = current; j < next; j += numChannels) {
                sum += fabs(m_buffer[absoluteOffset + j + channel]);
            }
        }
        squeezedBuffer[i] = sum * numChannels / (next - current);
    }
    return squeezedBuffer;
}

std::vector<float> Wave::getStretchedBuffer(uint32_t offset, uint32_t stretchedSampleCount, float stretchFactor,
                                            int channel) const
{
    std::vector<float> stretchedBuffer(stretchedSampleCount);
    uint32_t sampleCount = uint32_t(roundf(stretchedSampleCount / stretchFactor));
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint32_t absoluteOffset = numChannels * offset;
    for (uint32_t i = channel; i < sampleCount * numChannels; i += numChannels) {
        uint32_t current = uint32_t(roundf(i / numChannels * stretchFactor));
        uint32_t next = uint32_t(roundf((i + numChannels) / numChannels * stretchFactor));
        for (uint32_t j = current; j < next; ++j) {
            if (j < stretchedSampleCount) {
                stretchedBuffer[j] = m_buffer[absoluteOffset + i + channel];
            }
        }
    }
    return stretchedBuffer;
}

float Wave::getAbsPeak(int channel) const
{
    if (getNumFrames() > 1) {
        uint16_t numChannels = m_waveProperties.getNumChannels();
        float max = m_buffer[channel];
        for (uint32_t i = channel + numChannels; i < m_numSamples; i += numChannels) {
            if (fabs(m_buffer[i]) > max) {
                max = fabs(m_buffer[i]);
            }
        }
        return max;
    }
    return 0.f;
}

void Wave::insertAudio(uint32_t offset, const float* audio, uint32_t numSamples)
{
    detach();
    std::memcpy(m_buffer.get() + offset, audio, numSamples * sizeof(float));
}

void Wave::insertAudio(uint32_t offset, std::vector<float>& audio)
{
    insertAudio(offset, audio.data(), uint32_t(audio.size()));
}

bool Wave::isEmpty() const
{
    return (!m_buffer || m_numSamples == 0);
}

Wave& Wave::append(const Wave& other)
{
    if (other.getNumChannels() == getNumChannels()) {
        uint32_t numSamples = m_numSamples;
        resizeMemory(numSamples + other.getNumSamples());
        copySamples(other.m_buffer.get(), m_buffer.get(), other.m_numSamples, 0, numSamples);
    } else if (!isMono() && other.isMono()) {
        uint32_t numSamples = m_numSamples;
        uint16_t numChannels = getNumChannels();
        resizeMemory(numSamples + other.getNumSamples() * numChannels);
        for (uint32_t y = numSamples; y < m_numSamples; y += numChannels) {
            for (uint16_t x = 0; x < numChannels; ++x) {
                m_buffer[y + x] = other.m_buffer[(y - numSamples) / numChannels];
            }
        }
    }
    return *this;
}

float Wave::avgValue(int channel) const
{
    uint16_t numChannels = m_waveProperties.getNumChannels();
    double sum = 0.;
    for (uint32_t i = channel + numChannels; i < m_numSamples; i += numChannels) {
        sum += double(m_buffer[i]);
    }
    return float(sum / getNumFrames());
}

void Wave::changeVolume(float volume, int channel)
{
    detach();
    uint16_t numChannels = m_waveProperties.getNumChannels();
    float max = getAbsPeak(channel);
    float factor = volume / max;
    for (uint32_t i = channel; i < m_numSamples; i += numChannels) {
        m_buffer[i] *= factor;
    }
}

void Wave::downmixToMono()
{
    uint32_t numSamples = getNumSamples();
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint32_t newDataChunkSize = m_waveProperties.getDataChunkSize() / numChannels;
    uint32_t newNumSamples = m_numSamples / numChannels;
    auto newBuffer = std::shared_ptr<float[]>(new float[newNumSamples]);
    for (uint32_t y = 0; y < numSamples; y += numChannels) {
        float sum = 0.f;
        for (uint16_t x = 0; x < numChannels; ++x) {
            sum += m_buffer[y + x];
        }
        newBuffer[y / numChannels] = sum / float(numChannels);
    }
    m_buffer = std::move(newBuffer);
    m_waveProperties.setBlockAlign(m_waveProperties.getBlockAlign() / numChannels);
    m_waveProperties.setNumBytesPerSecond(m_waveProperties.getNumBytesPerSecond() / numChannels);
    m_waveProperties.setRiffChunkSize(newDataChunkSize + 36);
    m_waveProperties.setDataChunkSize(newDataChunkSize);
    m_waveProperties.setNumChannels(1);
    m_numSamples = newNumSamples;
}

Wave Wave::randomFromDuration(float duration, uint16_t bitDepth, uint32_t sampleRate)
{
    return generateRandom(uint32_t(roundf(sampleRate * duration)), bitDepth, sampleRate);
}

Wave Wave::generateRandom(uint32_t numFrames, uint16_t bitDepth, uint32_t sampleRate)
{
    Wave result(numFrames, 1, bitDepth, sampleRate);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1.f, 1.f - std::numeric_limits<float>::epsilon());
    for (uint32_t i = 0; i < numFrames; ++i) {
        result.m_buffer[i] = dis(gen);
    }
    return result;
}

Wave Wave::generateSine(float waveFreq, float phaseShift, uint32_t numFrames,
    uint16_t bitDepth, uint32_t sampleRate, bool multiThreaded)
{
    Wave result(numFrames, 1, bitDepth, sampleRate);
    constexpr float coeff = 1.f - std::numeric_limits<float>::epsilon();
    float timeStep = 1.f / sampleRate;
    float omega = 2 * 3.1415927f * waveFreq;
    #pragma omp parallel for if(multiThreaded)
    for (int32_t i = 0; i < int32_t(numFrames); ++i) {
        result.m_buffer[i] = coeff * sinf(omega * (i * timeStep + phaseShift));
    }
    return result;
}

Wave Wave::generateSquare(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames,
                          bool multiThreaded)
{
    Wave result(numFrames, 1, 16, samplingFreq);
    float timeStep = 1.f / samplingFreq;
    float omega = 2 * 3.1415927f * waveFreq;
    constexpr float coeff = 1.f - std::numeric_limits<float>::epsilon();
    #pragma omp parallel for if(multiThreaded)
    for (int32_t i = 0; i < int32_t(numFrames); ++i) {
        result.m_buffer[i] = sinf(omega * (i * timeStep + phaseShift)) < 0 ? -1.f * coeff : 1.f * coeff;
    }
    return result;
}

Wave Wave::generateTriangle(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames,
                            bool multiThreaded)
{
    Wave result(numFrames, 1, 16, samplingFreq);
    float timeStep = 1.f / samplingFreq;
    constexpr float pi = 3.1415927f;
    constexpr float coeff = 2 * (1.f - std::numeric_limits<float>::epsilon()) / pi;
    float omega = 2 * pi * waveFreq;
    #pragma omp parallel for if(multiThreaded)
    for (int32_t i = 0; i < int32_t(numFrames); ++i) {
        result.m_buffer[i] = coeff * asinf(sinf(omega * (i * timeStep + phaseShift)));
    }
    return result;
}

Wave Wave::generateClick(float bpm, uint8_t beatsPerBar, uint32_t numFrames,
                         uint16_t bitDepth, uint32_t sampleRate)
{
    Wave result(numFrames, 1, bitDepth, sampleRate);
    constexpr float pi = 3.1415927f;
    const double beatFrames = sampleRate * 60.0 / bpm;
    // Downbeat: 1000 Hz, 30 ms, amplitude 0.8; off-beat: 600 Hz, 20 ms, amplitude 0.5
    const uint32_t downLen = static_cast<uint32_t>(sampleRate * 0.030);
    const uint32_t beatLen = static_cast<uint32_t>(sampleRate * 0.020);

    double beatPos = 0.0;
    for (int beatNum = 0; beatPos < numFrames; ++beatNum) {
        const uint32_t start    = static_cast<uint32_t>(beatPos);
        const bool     isDown   = (beatNum % beatsPerBar == 0);
        const float    freq     = isDown ? 1000.f : 600.f;
        const float    amp      = isDown ? 0.8f   : 0.5f;
        const uint32_t clickLen = isDown ? downLen : beatLen;

        for (uint32_t i = 0; i < clickLen && start + i < numFrames; ++i) {
            const float t = static_cast<float>(i) / sampleRate;
            result.m_buffer[start + i] = amp * expf(-50.f * t) * sinf(2.f * pi * freq * t);
        }
        beatPos += beatFrames;
    }
    return result;
}

std::vector<float> Wave::getAveragedOutData(uint32_t binSize, bool absolute, int channel) const
{
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint32_t numFrames = getNumFrames();
    uint32_t numPoints = numFrames / binSize;
    uint32_t remainder = numFrames % binSize;
    if (remainder != 0) {
        ++numPoints;
    }
    std::vector<float> result(numPoints);
    for (uint32_t i = 0; i < numPoints; ++i) {
        double sum = 0.;
        uint32_t endPoint = 0;
        if (i == numPoints - 1) {
            endPoint = remainder == 0 ? binSize * numChannels : remainder * numChannels;
        } else {
            endPoint = binSize * numChannels;
        }
        if (!absolute) {
            for (uint32_t j = channel; j < endPoint; j += numChannels) {
                sum += double(m_buffer[i * binSize * numChannels + j]);
            }
        } else {
            for (uint32_t j = channel; j < endPoint; j += numChannels) {
                sum += double(std::fabs(m_buffer[i * binSize * numChannels + j]));
            }
        }
        result[i] = float(sum / binSize);
    }
    return result;
}

bool Wave::isLittleEndian() const
{
    return m_isLittleEndian;
}

bool Wave::isMono() const
{
    return (getNumChannels() == 1);
}

float Wave::maxValue(int channel) const
{
    uint16_t numChannels = m_waveProperties.getNumChannels();
    float max = 0.f;
    if (getNumFrames() > 1) {
        max = m_buffer[channel];
        for (uint32_t i = channel + numChannels; i < m_numSamples; i += numChannels) {
            if (m_buffer[i] > max) {
                max = m_buffer[i];
            }
        }
    }
    return max;
}

float Wave::minValue(int channel) const
{
    uint16_t numChannels = m_waveProperties.getNumChannels();
    float min = 0.f;
    if (getNumFrames() > 1) {
        min = m_buffer[channel];
        for (uint32_t i = channel + numChannels; i < m_numSamples; i += numChannels) {
            if (m_buffer[i] < min) {
                min = m_buffer[i];
            }
        }
    }
    return min;
}

void Wave::reverse(int channel)
{
    detach();
    dsp::reverseChannel(*this, channel);
}

Wave Wave::resample(uint32_t targetSampleRate) const
{
    if (targetSampleRate == 0) {
        throwError("Target sample rate must be greater than zero.",
                   "Wave::resample(uint32_t)");
    }
    if (isEmpty() || targetSampleRate == getSampleRate()) {
        return Wave(*this);
    }

    const uint32_t srcFrames = getNumFrames();
    const uint16_t channels  = getNumChannels();
    const uint32_t dstFrames = std::max<uint32_t>(
        1,
        static_cast<uint32_t>(std::llround(
            static_cast<double>(srcFrames) * targetSampleRate / getSampleRate())));

    Wave result(dstFrames, channels, getSampleBitDepth(), targetSampleRate);
    float* dst = result.audioData();

    if (srcFrames == 1) {
        for (uint32_t i = 0; i < dstFrames; ++i) {
            for (uint16_t ch = 0; ch < channels; ++ch) {
                dst[i * channels + ch] = m_buffer[ch];
            }
        }
        return result;
    }

    const double srcRate = static_cast<double>(getSampleRate());
    const double dstRate = static_cast<double>(targetSampleRate);
    for (uint32_t i = 0; i < dstFrames; ++i) {
        const double srcPos = i * srcRate / dstRate;
        const uint32_t i0 = static_cast<uint32_t>(srcPos);
        const uint32_t i1 = std::min(i0 + 1, srcFrames - 1);
        const float frac = static_cast<float>(srcPos - i0);

        for (uint16_t ch = 0; ch < channels; ++ch) {
            const float a = m_buffer[i0 * channels + ch];
            const float b = m_buffer[i1 * channels + ch];
            dst[i * channels + ch] = a + (b - a) * frac;
        }
    }

    return result;
}

void Wave::setAudio(const float* audio, uint32_t numSamples)
{
    if (numSamples != m_numSamples) {
        resizeMemory(numSamples, false);
    } else {
        detach();
    }
    if (numSamples == 0 || audio == nullptr) {
        return;
    }
    std::memcpy(m_buffer.get(), audio, m_numSamples * sizeof(float));
}

void Wave::setAudio(std::vector<float>& audio)
{
    setAudio(audio.data(), uint32_t(audio.size()));
}

void Wave::setLittleEndian(bool isLittleEndian)
{
    m_isLittleEndian = isLittleEndian;
}

void Wave::swapChannels(int from, int to)
{
    detach();
    uint32_t numSamples = getNumSamples();
    uint16_t numChannels = getNumChannels();
    if (numChannels < 2) {
        throwError("The number of channels must be at least two.",
                   "Wave::swapChannels()");
    }
    int idx1 = std::min(from, to);
    int idx2 = std::max(from, to);
    int offset = idx2 - idx1;
    for (uint32_t i = idx1; i < numSamples; i += numChannels) {
        float sample = m_buffer[i];
        m_buffer[i] = m_buffer[i + offset];
        m_buffer[i + offset] = sample;
    }
}

void Wave::upmixToStereo()
{
    uint32_t numSamples = getNumSamples();
    uint16_t numChannels = 2;
    uint32_t newDataChunkSize = m_waveProperties.getDataChunkSize() * numChannels;
    uint32_t newNumSamples = m_numSamples * numChannels;
    std::cout << newNumSamples << std::endl;
    auto newBuffer = std::shared_ptr<float[]>(new float[newNumSamples]);
    for (uint32_t y = 0; y < numSamples; ++y) {
        for (uint16_t x = 0; x < numChannels; ++x) {
            newBuffer[y * numChannels + x] = m_buffer[y];
        }
    }
    m_buffer = std::move(newBuffer);
    m_waveProperties.setBlockAlign(m_waveProperties.getBlockAlign() * numChannels);
    m_waveProperties.setNumBytesPerSecond(m_waveProperties.getNumBytesPerSecond() * numChannels);
    m_waveProperties.setRiffChunkSize(newDataChunkSize + 36);
    m_waveProperties.setDataChunkSize(newDataChunkSize);
    m_waveProperties.setNumChannels(numChannels);
    m_numSamples = newNumSamples;
}

void Wave::zeroInitHeader()
{
    std::memset(&m_header, 0, sizeof(Header));
    std::memset(&m_dataSubChunk, 0, sizeof(DataSubChunk));
}

uint32_t Wave::getDataChunkSize() const
{
    return m_waveProperties.getDataChunkSize();
}

uint16_t Wave::getNumChannels() const
{
    return m_waveProperties.getNumChannels();
}

uint32_t Wave::getNumFrames() const
{
    return m_numSamples / m_waveProperties.getNumChannels();
}

uint32_t Wave::getNumSamples() const
{
    return m_numSamples;
}

uint16_t Wave::getSampleBitDepth() const
{
    return m_waveProperties.getNumBitsPerSample();
}

uint32_t Wave::getSampleRate() const
{
    return m_waveProperties.getSamplingFrequency();
}

void Wave::saveAs(const char* filename, uint32_t bufferSize)
{
    FILE* pFile = std::fopen(filename, "wb");
    if (!pFile) {
        throwError("Failed to create the file.",
                   "Wave::saveAs(const char*, uint32_t)");
    }
    generateHeader();
    std::fwrite(&m_header, 1, sizeof(Header), pFile);
    std::fwrite(&m_dataSubChunk, 1, sizeof(DataSubChunk), pFile);
    writeData(pFile, bufferSize);
    std::fclose(pFile);
}

void Wave::saveAs(const std::string& filename, uint32_t bufferSize)
{
    saveAs(filename.c_str(), bufferSize);
}

void Wave::writeData(std::FILE* file, uint32_t bufferSize)
{
    uint8_t* pBuffer = reinterpret_cast<uint8_t*>(std::malloc(bufferSize));
    if (pBuffer == nullptr) {
        throwError("Failed to create a buffer.",
                   "Wave::writeData(std::FILE*, uint32_t)");
    }

    uint16_t sampleBitDepth = m_waveProperties.getNumBitsPerSample();
    uint32_t dataChunkSize = m_waveProperties.getDataChunkSize();
    uint32_t numSamples = dataChunkSize / (uint32_t(sampleBitDepth) / 8);
    uint16_t numChannels = m_waveProperties.getNumChannels();
    uint16_t frameSize = numChannels * sampleBitDepth / 8;
    if (bufferSize % frameSize != 0) {
        throwError("The buffer size must be a multiple of the frame size.",
                   "Wave::readData(std::FILE*, uint32_t)");
    }

    const float* src = m_buffer.get();
    bool isEndiannessMismatched = isCpuBigEndian() == m_isLittleEndian;
    uint32_t offset = 0;
    if (sampleBitDepth == 8) {
        while (offset < numSamples) {
            uint32_t numBytesToWrite = std::min(bufferSize, dataChunkSize - offset);
            uint32_t numBufferFrames = numBytesToWrite / frameSize;
            for (uint32_t y = 0; y < numBufferFrames; ++y) {
                for (uint16_t x = 0; x < numChannels; ++x) {
                    int i = y * numChannels + x;
                    pBuffer[i] = uint8_t(roundf((src[offset + i] + 1.f) * 127.5f));
                }
            }
            std::fwrite(pBuffer, 1, numBytesToWrite, file);
            offset += bufferSize;
        }
    } else if (sampleBitDepth == 16) {
        while (offset < numSamples) {
            uint32_t numBytesToWrite = std::min(bufferSize, dataChunkSize - 2 * offset);
            uint32_t numBufferFrames = numBytesToWrite / frameSize;
            int16_t* pI16Buffer = reinterpret_cast<int16_t*>(pBuffer);
            for (uint32_t y = 0; y < numBufferFrames; ++y) {
                for (uint16_t x = 0; x < numChannels; ++x) {
                    uint32_t i = y * numChannels + x;
                    pI16Buffer[i] = int16_t(src[offset + i] * 32768.f);
                }
            }
            if (isEndiannessMismatched) {
                changeBufferEndianness(pBuffer, 2, numBytesToWrite);
            }
            std::fwrite(pBuffer, 1, numBytesToWrite, file);
            offset += bufferSize / 2;
        }
    } else if (sampleBitDepth == 24) {
        while (offset < numSamples) {
            uint32_t numBytesToWrite = std::min(bufferSize, dataChunkSize - 3 * offset);
            uint32_t numBufferFrames = numBytesToWrite / frameSize;
            for (uint32_t y = 0; y < numBufferFrames; ++y) {
                for (uint16_t x = 0; x < numChannels; ++x) {
                    uint32_t i = 3 * (y * numChannels + x);
                    int val = int(src[offset + i / 3] * 2147483648.f);
                    pBuffer[i] = uint8_t(val >> 8);
                    pBuffer[i + 1] = uint8_t(val >> 16);
                    pBuffer[i + 2] = uint8_t(val >> 24);
                }
            }
            if (isEndiannessMismatched) {
                changeBufferEndianness(pBuffer, 3, numBytesToWrite);
            }
            std::fwrite(pBuffer, 1, numBytesToWrite, file);
            offset += bufferSize / 3;
        }
    } else if (sampleBitDepth == 32) {
        while (offset < numSamples) {
            uint32_t numBytesToWrite = std::min(bufferSize, dataChunkSize - 4 * offset);
            uint32_t numBufferFrames = numBytesToWrite / frameSize;
            float* pF32Buffer = reinterpret_cast<float*>(pBuffer);
            copySamples(src, pF32Buffer, numBufferFrames * numChannels, offset, 0);
            if (isEndiannessMismatched) {
                changeBufferEndianness(pBuffer, 4, numBytesToWrite);
            }
            std::fwrite(pBuffer, 1, numBytesToWrite, file);
            offset += bufferSize / 4;
        }
    } else {
        std::free(pBuffer);
        throwError("Bit depth not supported.",
                   "Wave::writeData(std::FILE*, uint32_t)");
    }
    std::free(pBuffer);
}

void Wave::setSampleBitDepth(uint16_t sampleBitDepth)
{
    if (sampleBitDepth != 8 && sampleBitDepth != 16 && sampleBitDepth != 24 && sampleBitDepth != 32) {
        throwError("Bit depth not supported.",
                   "Wave::setSampleBitDepth(uint32_t)");
    }
    m_waveProperties.setNumBitsPerSample(sampleBitDepth);
    m_waveProperties.setAudioFormat(sampleBitDepth != 32 ? 1 : 3);
    uint16_t numChannels = m_waveProperties.getNumChannels();
    m_waveProperties.setNumBytesPerSecond(numChannels * m_waveProperties.getSamplingFrequency() * sampleBitDepth / 8);
    m_waveProperties.setBlockAlign(numChannels * sampleBitDepth / 8);
    m_waveProperties.setDataChunkSize(m_numSamples * sampleBitDepth / 8);
}

void Wave::setSampleRate(uint32_t sampleRate)
{
    m_waveProperties.setSamplingFrequency(sampleRate);
}

float& Wave::operator()(uint32_t sample, int channel)
{
    detach();
    return m_buffer[sample * m_waveProperties.getNumChannels() + channel];
}

const float& Wave::operator()(uint32_t sample, int channel) const
{
    return m_buffer[sample * m_waveProperties.getNumChannels() + channel];
}

void Wave::operator=(const Wave& other)
{
    if (this == &other) {
        return;
    }
    m_header        = other.m_header;
    m_dataSubChunk  = other.m_dataSubChunk;
    m_isLittleEndian = other.m_isLittleEndian;
    m_numSamples    = other.m_numSamples;
    m_waveProperties = other.m_waveProperties;
    m_buffer        = other.m_buffer;
}

std::ostream& operator<<(std::ostream& os, const Wave& wav)
{
    os << "ChunkID: " << std::string(wav.m_header.riff.descriptor.id)
        .substr(0, sizeof(wav.m_header.riff.descriptor.id)) << std::endl;
    os << "ChunkSize: " << wav.m_waveProperties.getRiffChunkSize() << std::endl;
    os << "Format: " << std::string(wav.m_header.riff.type)
        .substr(0, sizeof(wav.m_header.riff.type)) << std::endl;
    os << "----------" << std::endl;
    os << "Subchunk1ID: " << std::string(wav.m_header.wave.descriptor.id)
        .substr(0, sizeof(wav.m_header.wave.descriptor.id)) << std::endl;
    os << "Subchunk1Size: " << wav.m_waveProperties.getFmtChunkSize() << std::endl;
    os << "AudioFormat: " << wav.m_waveProperties.getAudioFormat() << std::endl;
    os << "NumChannels: " << wav.m_waveProperties.getNumChannels() << std::endl;
    os << "SampleRate: " << wav.m_waveProperties.getSamplingFrequency() << std::endl;
    os << "ByteRate: " << wav.m_waveProperties.getNumBytesPerSecond() << std::endl;
    os << "BlockAlign: " << wav.m_waveProperties.getBlockAlign() << std::endl;
    os << "BitsPerSample: " << wav.m_waveProperties.getNumBitsPerSample() << std::endl;
    os << "----------" << std::endl;
    os << "Subchunk2ID: " << std::string(wav.m_dataSubChunk.descriptor.id)
        .substr(0, sizeof(wav.m_dataSubChunk.descriptor.id)) << std::endl;
    os << "Subchunk2Size: " << wav.m_waveProperties.getDataChunkSize() << std::endl << std::endl;
    return os;
}

template<typename T>
inline T Wave::reverseBytes(T val)
{
    char* pVal = reinterpret_cast<char*>(&val);
    T res = 0;
    char* pRes = reinterpret_cast<char*>(&res);
    size_t count = sizeof(T);
    std::reverse_copy(pVal, pVal + count, pRes);
    return res;
}
}
