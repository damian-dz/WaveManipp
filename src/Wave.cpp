#include "Wave.hpp"
#include <iostream>

namespace wm {

    /*!
     * \brief Constructs an empty Wave object with a zero-initialized header.
     */
    Wave::Wave():
        m_pData(nullptr),
        m_dataSize(0)
    {
        zeroInitHeader();
    }

    Wave::Wave(const char *c_pFilename)
    {
        zeroInitHeader();
        open(c_pFilename);
    }

    Wave::Wave(const std::string &filename)
    {
        zeroInitHeader();
        open(filename);
    }

    /*!
     * \brief Destroys the Wave object and frees its associated resources.
     */
    Wave::~Wave()
    {
        if (m_pData != nullptr) {
            std::free(m_pData);
            m_pData = nullptr;
        }
    }

    void Wave::generateHeader(uint32_t dataSize, uint16_t numChannels, uint32_t sampleRate, uint16_t bitsPerSample)
    {
        if (bitsPerSample != 8 && bitsPerSample != 16 && bitsPerSample != 24 && bitsPerSample != 32) {
            throwError("Bit depth not supported.",
                       "Wave::generateHeader(uint32_t, uint16_t, uint32_t, uint16_t)");
        }

        std::memcpy(m_combinedHeader.riff.descriptor.id, "RIFF", 4);
        m_combinedHeader.riff.descriptor.size = dataSize + 36;
        std::memcpy(m_combinedHeader.riff.type, "WAVE", 4);

        std::memcpy(m_combinedHeader.wave.descriptor.id, "fmt ", 4);
        m_combinedHeader.wave.descriptor.size = 16;
        m_combinedHeader.wave.audioFormat = bitsPerSample != 32 ? 1 : 3;
        m_combinedHeader.wave.numChannels = numChannels;
        m_combinedHeader.wave.sampleRate = sampleRate;
        m_combinedHeader.wave.byteRate = sampleRate * numChannels * bitsPerSample / 8;
        m_combinedHeader.wave.blockAlign = numChannels * bitsPerSample / 8;
        m_combinedHeader.wave.bitsPerSample = bitsPerSample;

        std::memcpy(m_dataHeader.descriptor.id, "data", 4);
        m_dataHeader.descriptor.size = dataSize;
    }

    void setLimits(int& x0, int& xLim, Channels channels)
    {
        if (channels == Both) {
            x0 = 0;
        } else if (channels == Left) {
            x0 = 0;
            xLim = 1;
        } else if (channels == Right) {
            x0 = 1;
            xLim = 2;
        }
    }

    bool peekForId(const std::string &id, std::FILE *pFile)
    {
        bool res = true;
        int pos = std::ftell(pFile);
        for (unsigned i = 0; i < id.size(); ++i) {
            int c = std::fgetc(pFile);
            if (id[i] != (char)c) {
                res = false;
                break;
            }
        }
        std::fseek(pFile, pos, 0);
        return res;
    }

    void findDataChunk(std::FILE *pFile)
    {
        int pos = std::ftell(pFile) + 2;
        std::fseek(pFile, pos, 0);
        std::string idData = "data";
        while (!peekForId(idData, pFile)) {
            pos += 2;
            std::fseek(pFile, pos, 0);
        }
    }

    /*!
     * \brief Loads the file specified by its name into memory.
     * \param c_pFilename &mdash; the name of the file
     * \param bufferSize &mdash; the size of the intermediate buffer (default 24576)
     */
    void Wave::open(const char *c_pFilename, size_t bufferSize)
    {
        std::FILE *pFile = std::fopen(c_pFilename, "rb");
        if (!pFile) {
            throwError("Failed to open the file.",
                "Wave::open(const std::string &, size_t)");
        }
        std::fread(&m_combinedHeader, 1, sizeof(CombinedHeader), pFile);

        if (m_combinedHeader.wave.descriptor.size > 16) {
            int pos = std::ftell(pFile);
            std::fseek(pFile, pos + (m_combinedHeader.wave.descriptor.size - 16), 0);
        }
        if (!peekForId(std::string("data"), pFile)) {
            findDataChunk(pFile);
        }
        std::fread(&m_dataHeader, 1, sizeof(DATAHeader), pFile);
        m_dataSize = m_dataHeader.descriptor.size;
        m_numSamples = m_dataSize / (int(m_combinedHeader.wave.bitsPerSample) / 8);
        m_pData = reinterpret_cast<float*>(std::malloc(m_dataSize * sizeof(m_pData)));
        if (m_pData != nullptr) {
            readData(pFile, bufferSize);
        }
        std::fclose(pFile);
    }

    /*!
     * \brief Loads the file specified by its name into memory.
     * \param filename &mdash; the name of the file
     * \param bufferSize &mdash; the size of the intermediate buffer (default 24576)
     */
    void Wave::open(const std::string &filename, size_t bufferSize)
    {
        open(filename.c_str(), bufferSize);
    }

    void Wave::readData(std::FILE *pFile, size_t bufferSize)
    {
        uint8_t *pBuffer = reinterpret_cast<uint8_t *>(std::malloc(bufferSize));
        if (pBuffer == nullptr) {
            throwError("Failed to create a buffer.",
                       "Wave::readData(std::FILE *, size_t)");
        }

        uint16_t sampleBitDepth = getSampleBitDepth();
        uint32_t dataSize = getDataSize();
        int numSamples = dataSize / (int(sampleBitDepth) / 8);
        int numChannels = getNumChannels();
        int frameSize = numChannels * sampleBitDepth / 8;
        if (bufferSize % frameSize != 0) {
            throwError("The buffer size must be a multiple of the frame size.",
                       "Wave::readData(std::FILE *, size_t)");
        }

        int numFrames = numSamples / numChannels;
        size_t offset = 0;
        if (sampleBitDepth == 8) {
            while (offset < numSamples) {
                size_t numBytesToRead = std::min(bufferSize, dataSize - offset);
                std::fread(pBuffer, 1, numBytesToRead, pFile);
                int numBufferFrames = numBytesToRead / frameSize;
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        m_pData[offset + i] = static_cast<float>(pBuffer[i] / 127.5f - 1.f);
                    }
                }
                offset += bufferSize;
            }
        } else if (sampleBitDepth == 16) {
            while (offset < numSamples) {
                size_t numBytesToRead = std::min(bufferSize, dataSize - offset * 2);
                int numBufferFrames = numBytesToRead / frameSize;
                int16_t *pI16Buffer = reinterpret_cast<int16_t *>(pBuffer);
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        m_pData[offset + i] = static_cast<float>(pI16Buffer[i] / 32768.f);
                    }
                }
                offset += bufferSize / 2;
            }
        } else if (sampleBitDepth == 24) {
            while (offset < numSamples) {
                size_t numBytesToRead = std::min(bufferSize, dataSize - offset * 3);
                std::fread(pBuffer, 1, numBytesToRead, pFile);
                int numBufferFrames = numBytesToRead / frameSize;
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = (y * numChannels + x) * 3;
                        m_pData[offset + i / 3] = static_cast<float>(
                            (pBuffer[i] << 8 | pBuffer[i + 1] << 16 | pBuffer[i + 2] << 24) / 2147483648.f);
                    }
                }
                offset += bufferSize / 3;
            }     
        } else if (sampleBitDepth == 32) {
            while (offset < numSamples) {
                size_t numBytesToRead = std::min(bufferSize, dataSize - offset * 4);
                std::fread(pBuffer, 1, numBytesToRead, pFile);
                int numBufferFrames = numBytesToRead / frameSize;
                float *pF32Buffer = reinterpret_cast<float *>(pBuffer);
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        m_pData[offset + i] = static_cast<float>(pF32Buffer[i]);
                    }
                }
                offset += bufferSize / 4;
            }
        } else {
            std::free(pBuffer);
            throwError("Bit depth not supported.",
                       "Wave::readData(std::FILE *, size_t)");
        }
        std::free(pBuffer);
    }

    std::vector<float> Wave::getBuffer(size_t offset, size_t sampleCount, Channels channels) const
    {
        std::vector<float> buffer(sampleCount);
        for (int i = 0; i < sampleCount; ++i) {
            buffer[i] = m_pData[offset + i];
        }
        return buffer;
    }

    void Wave::changeVolume(float volume, Channels channels)
    { 
        int numSamples = getNumSamples();
        int numChannels = getNumChannels();
        int numFrames = numSamples / numChannels;
        int x0 = 0;
        int xLim = numChannels;
        setLimits(x0, xLim, channels);
        float max[2] = { std::numeric_limits<float>::min(), std::numeric_limits<float>::min() };
        for (int y = 0; y < numFrames; ++y) {
            for (int x = x0; x < xLim; ++x) {
                if (fabs(m_pData[y * numChannels + x]) > max[x]) {
                    max[x] = fabs(m_pData[y * numChannels + x]);
                }
            }
        }
        float factors[2] = { volume / max[0], volume / max[1] };
        for (int y = 0; y < numFrames; ++y) {
            for (int x = x0; x < xLim; ++x) {
                m_pData[y * numChannels + x] *= factors[x];
            }
        }
    }

    void Wave::downmixToMono()
    {
        uint32_t numSamples = getNumSamples();
        uint16_t numChannels = getNumChannels();
        uint32_t newDataSize = m_dataSize / numChannels;
        uint32_t newNumSamples = m_numSamples / numChannels;
        float *pNewData = reinterpret_cast<float *>(std::malloc(newNumSamples * sizeof(pNewData)));
        if (pNewData != nullptr) {
            for (int i = 0; i < numSamples; i += numChannels) {
                pNewData[i >> 1] = (m_pData[i] + m_pData[i + 1]) / 2.f;
            }
            std::free(m_pData);
            m_pData = pNewData;
            m_combinedHeader.wave.blockAlign /= numChannels;
            m_combinedHeader.wave.byteRate /= numChannels;
            m_combinedHeader.riff.descriptor.size = newDataSize + 36;
            m_dataHeader.descriptor.size = newDataSize;
            m_dataSize = newDataSize;
            m_numSamples = newNumSamples;
            m_combinedHeader.wave.numChannels = 1;
        }
    }

    void Wave::swapChannels()
    {
        int numSamples = getNumSamples();
        if (getNumChannels() != 2) {
            throwError("The number of channels must be two.",
                "Wave::swapChannels()");
        }

        for (int i = 0; i < numSamples; i += 2) {
            float val = m_pData[i];
            m_pData[i] = m_pData[i + 1];
            m_pData[i + 1] = val;
        }
    }

    void Wave::zeroInitHeader()
    {
        std::memset(&m_combinedHeader, 0, sizeof(CombinedHeader));
        std::memset(&m_dataHeader, 0, sizeof(DATAHeader));
    }

    /*!
     * \brief Gets the number of bytes necessary to store the audio data with the current settings.
     * \result The data size.
     */
    uint32_t Wave::getDataSize() const
    {
        return m_dataHeader.descriptor.size;
    }

    /*!
     * \brief Gets the number of channels for the audio data.
     * \result The number of channels.
     */
    uint16_t Wave::getNumChannels() const
    {
        return m_combinedHeader.wave.numChannels;
    }

    uint32_t Wave::getNumSamples() const
    {
        return m_numSamples;
    }

    /*!
     * \brief Gets the number of bits currently used to represent a single audio sample.
     * \result The sample bit depth.
     */
    uint16_t Wave::getSampleBitDepth() const
    {
        return m_combinedHeader.wave.bitsPerSample;
    }

    void Wave::saveAs(const char *c_pFilename, uint32_t sampleRate, uint16_t sampleBitDepth, size_t bufferSize)
    {
        FILE *pFile = std::fopen(c_pFilename, "wb");
        if (!pFile) {
            throwError("Failed to create the file.",
                "Wave::saveAs(const char *, uint32_t, uint16_t, size_t)");
        }
        generateHeader(m_numSamples * sampleBitDepth / 8, getNumChannels(), sampleRate, sampleBitDepth);
        std::fwrite(&m_combinedHeader, 1, sizeof(CombinedHeader), pFile);
        std::fwrite(&m_dataHeader, 1, sizeof(DATAHeader), pFile);
        writeData(pFile, bufferSize);
        std::fclose(pFile);
    }

    void Wave::saveAs(const std::string &filename, uint32_t sampleRate, uint16_t sampleBitDepth, size_t bufferSize)
    {
        saveAs(filename.c_str(), sampleRate, sampleBitDepth, bufferSize);
    }

    void Wave::writeData(std::FILE *pFile, size_t bufferSize)
    {
        uint8_t *pBuffer = reinterpret_cast<uint8_t *>(std::malloc(bufferSize));
        if (pBuffer == nullptr) {
            throwError("Failed to create a buffer.",
                       "Wave::writeData(std::FILE *, size_t)");
        }

        uint16_t sampleBitDepth = this->getSampleBitDepth();
        uint32_t dataSize = this->getDataSize();
        int numSamples = dataSize / (int(sampleBitDepth) / 8);
        int numChannels = this->getNumChannels();
        int frameSize = numChannels * sampleBitDepth / 8;
        if (bufferSize % frameSize != 0) {
            throwError("The buffer size must be a multiple of the frame size.",
                "Wave::readData(std::FILE *, size_t)");
        }

        int numFrames = numSamples / numChannels;
        size_t offset = 0;
        if (sampleBitDepth == 8) {
            while (offset < numSamples) {
                size_t numBytesToWrite = std::min(bufferSize, dataSize - offset);
                int numBufferFrames = numBytesToWrite / frameSize;
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        pBuffer[i] = uint8_t(roundf((m_pData[offset + i] + 1.f) * 127.5f));
                    }
                }
                std::fwrite(pBuffer, 1, numBytesToWrite, pFile);
                offset += bufferSize;
            }
        } else if (sampleBitDepth == 16) {
            while (offset < numSamples) {
                size_t numBytesToWrite = std::min(bufferSize, dataSize - offset * 2);
                int numBufferFrames = numBytesToWrite / frameSize;
                int16_t *pI16Buffer = reinterpret_cast<int16_t *>(pBuffer);
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        pI16Buffer[i] = int16_t(m_pData[offset + i] * 32768.f);
                    }
                }
                std::fwrite(pBuffer, 1, numBytesToWrite, pFile);
                offset += bufferSize / 2;
            }
        } else if (sampleBitDepth == 24) {
            while (offset < numSamples) {
                size_t numBytesToWrite = std::min(bufferSize, dataSize - offset * 3);
                int numBufferFrames = numBytesToWrite / frameSize;
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = (y * numChannels + x) * 3;
                        int val = int(m_pData[offset + i / 3] * 2147483648.f);
                        pBuffer[i] = uint8_t(val >> 8);
                        pBuffer[i + 1] = uint8_t(val >> 16);
                        pBuffer[i + 2] = uint8_t(val >> 24);
                    }
                }
                std::fwrite(pBuffer, 1, numBytesToWrite, pFile);
                offset += bufferSize / 3;
            }
        } else if (sampleBitDepth == 32) {
            while (offset < numSamples) {
                size_t numBytesToWrite = std::min(bufferSize, dataSize - offset * 4);
                int numBufferFrames = numBytesToWrite / frameSize;
                float *pF32Buffer = reinterpret_cast<float *>(pBuffer);
                for (int y = 0; y < numBufferFrames; ++y) {
                    for (int x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        pF32Buffer[i] = m_pData[offset + i];
                    }
                }
                std::fwrite(pBuffer, 1, numBytesToWrite, pFile);
                offset += bufferSize / 4;
            }
        } else {
            std::free(pBuffer);
            throwError("Bit depth not supported.",
                       "Wave::writeData(std::FILE *, size_t)");
        }
        std::free(pBuffer);
    }

    std::ostream &operator<<(std::ostream &os, const Wave &wav)
    {
        os << "ChunkID: " << std::string(wav.m_combinedHeader.riff.descriptor.id)
            .substr(0, sizeof(wav.m_combinedHeader.riff.descriptor.id)) << std::endl;
        os << "ChunkSize: " << wav.m_combinedHeader.riff.descriptor.size << std::endl;
        os << "Format: " << std::string(wav.m_combinedHeader.riff.type)
            .substr(0, sizeof(wav.m_combinedHeader.riff.type)) << std::endl;
        os << "----------" << std::endl;
        os << "Subchunk1ID: " << std::string(wav.m_combinedHeader.riff.descriptor.id)
            .substr(0, sizeof(wav.m_combinedHeader.wave.descriptor.id)) << std::endl;
        os << "Subchunk1Size: " << wav.m_combinedHeader.wave.descriptor.size << std::endl;
        os << "AudioFormat: " << wav.m_combinedHeader.wave.audioFormat << std::endl;
        os << "NumChannels: " << wav.m_combinedHeader.wave.numChannels << std::endl;
        os << "SampleRate: " << wav.m_combinedHeader.wave.sampleRate << std::endl;
        os << "ByteRate: " << wav.m_combinedHeader.wave.byteRate << std::endl;
        os << "BlockAlign: " << wav.m_combinedHeader.wave.blockAlign << std::endl;
        os << "BitsPerSample: " << wav.m_combinedHeader.wave.bitsPerSample << std::endl;
        os << "----------" << std::endl;
        os << "Subchunk2ID: " << std::string(wav.m_dataHeader.descriptor.id)
            .substr(0, sizeof(wav.m_dataHeader.descriptor.id)) << std::endl;
        os << "Subchunk2Size: " << wav.m_dataHeader.descriptor.size << std::endl << std::endl;
        return os;
    }

}