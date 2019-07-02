#include "Wave.hpp"
#include <iostream>

namespace wm {

    /*!
     * \brief Constructs an empty Wave object with a zero-initialized header.
     */
    Wave::Wave():
        m_pData(nullptr)
    {
        zeroInitHeader();
    }

    /*!
     * \brief Loads the specified file into a Wave object in its entirety.
     * \param filename &mdash; the name of the file
     */
    Wave::Wave(const char* filename)
    {
        zeroInitHeader();
        open(filename);
    }

    /*!
     * \brief Loads the specified file into a Wave object in its entirety.
     * \param filename &mdash; the name of the file
     */
    Wave::Wave(const std::string& filename)
    {
        zeroInitHeader();
        open(filename);
    }

    Wave::Wave(uint32_t numFrames, uint16_t numChannels, uint16_t bitDepth, uint32_t sampleRate)
    {
        zeroInitHeader();
        setFourCharacterCodes();
        m_waveProperties.setNumChannels(numChannels);
        m_isLittleEndian = true;
        m_numSamples = numChannels * numFrames;
        reserveMemory(m_numSamples);
        setSampleRate(sampleRate);
        setSampleBitDepth(bitDepth);
        m_waveProperties.setRiffChunkSize(m_waveProperties.getDataChunkSize() + 36);
        m_waveProperties.setFmtChunkSize(16);
    }

    /*!
     * \brief Creates a deep copy of the input Wave object.
     * \param other &mdash; the object to copy
     */
    Wave::Wave(const Wave& other) :
        m_header(other.m_header),
        m_dataSubChunk(other.m_dataSubChunk),
        m_waveProperties(other.m_waveProperties),
        m_numSamples(other.m_numSamples),
        m_isLittleEndian(other.m_isLittleEndian)
    {
        if (other.m_pData != nullptr) {
            reserveMemory(other.m_numSamples);
            copySamples(other.m_pData, m_pData, m_numSamples);
        } else {
            m_pData = nullptr;
        }    
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

    void Wave::changeBufferEndianness(uint8_t* buffer, uint32_t sampleLength, uint32_t bufferLength)
    {
        for (uint32_t i = 0; i < bufferLength; i += sampleLength) {
            uint32_t j = i + sampleLength - 1;
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

    /*!
     * \brief Reserves enough memory to store the specified number of samples using floating-point format.
     * \param numSamples &mdash; the number of samples
     * \param zeroInit &mdash; specifies if memory should be set to zero
     */
    void Wave::reserveMemory(uint32_t numSamples, bool zeroInit)
    {
        if (!zeroInit) {
            m_pData = reinterpret_cast<float*>(std::malloc(numSamples * sizeof(m_pData)));
        } else {
            m_pData = reinterpret_cast<float*>(std::calloc(numSamples, sizeof(m_pData)));
        }
    }

    /*!
     * \brief Resizes the block of memory to accomodate the specified number of samples using floating-point format.
     * \param numSamples &mdash; the new number of samples
     * \param zeroInit &mdash; specifies if any additional memory should be set to zero
     */
    void Wave::resizeMemory(uint32_t numSamples, bool zeroInit)
    {
        if (numSamples == m_numSamples) {
            throwError("The number of samples remains the same. Resize unnecessary.",
                "Wave::resizeMemory(uint32_t, bool)");
        }
        uint32_t count = std::min(m_numSamples, numSamples);
        m_numSamples = numSamples;
        m_waveProperties.setDataChunkSize(numSamples * m_waveProperties.getNumBitsPerSample() / 8);
        if (numSamples == 0) {
            std::free(m_pData);
            m_pData = nullptr;
        } else {
            float* pNewData = reinterpret_cast<float*>(std::realloc(m_pData, numSamples * sizeof(m_pData)));
            if (pNewData == nullptr) {
                throwError("Failed to reallocate memory.",
                    "Wave::resizeMemory(uint32_t, bool)");
            }
            m_pData = pNewData;
            if (zeroInit) {
                for (uint32_t i = count; i < numSamples; ++i) {
                    m_pData[i] = 0.f;
                }
            }
        }
    }

    void Wave::copySamples(float* source, float* destination, uint32_t count, uint32_t srcOffset, uint32_t destOffset)
    {
        std::memcpy(&destination[destOffset], &source[srcOffset], count * sizeof(float));
    }

    bool Wave::peekForId(const std::string &id, std::FILE *pFile)
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

    void Wave::findDataChunk(std::FILE *pFile)
    {
        int pos = std::ftell(pFile) + 2;
        std::fseek(pFile, pos, 0);
        std::string id = "data";
        while (!peekForId(id, pFile)) {
            pos += 2;
            std::fseek(pFile, pos, 0);
        }
    }

    /*!
     * \brief Loads the file specified by its name into memory.
     * \param filename &mdash; the name of the file
     * \param bufferSize &mdash; the size of the intermediate read buffer (default 24576)
     */
    void Wave::open(const char* filename, uint32_t bufferSize)
    {
        std::FILE* pFile = std::fopen(filename, "rb");
        if (!pFile) {
            throwError("Failed to open the file.",
                "Wave::open(const std::string&, uint32_t)");
        }
        std::fread(&m_header, 1, sizeof(Header), pFile);

        if (std::memcmp(m_header.riff.descriptor.id, "RIFF", 4) == 0) {
            m_isLittleEndian = true;
        } else if (std::memcmp(m_header.riff.descriptor.id, "RIFX", 4) == 0) {
            m_isLittleEndian = false;
        } else {
            throwError("Wrong format or file corrupt.",
                "Wave::open(const std::string&, uint32_t)");
        }
        
        setWaveProperties();

        if (m_waveProperties.getFmtChunkSize() > 16) {
            int pos = std::ftell(pFile);
            std::fseek(pFile, pos + (m_waveProperties.getFmtChunkSize() - 16), 0);
        }
        if (!peekForId(std::string("data"), pFile)) {
            findDataChunk(pFile);
        }
        std::fread(&m_dataSubChunk, 1, sizeof(DataSubChunk), pFile);

        setWaveProperties(true);

        m_numSamples = m_waveProperties.getDataChunkSize() / (m_waveProperties.getNumBitsPerSample() / 8);
        reserveMemory(m_numSamples);
        if (m_pData != nullptr) {
            readData(pFile, bufferSize);
        }
        std::fclose(pFile);
    }

    /*!
     * \brief Loads the file specified by its name into memory.
     * \param filename &mdash; the name of the file
     * \param bufferSize &mdash; the size of the intermediate read buffer (default 24576)
     */
    void Wave::open(const std::string &filename, uint32_t bufferSize)
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
        uint32_t numSamples = dataChunkSize / (uint32_t(sampleBitDepth) / 8);
        uint16_t numChannels = m_waveProperties.getNumChannels();
        uint16_t frameSize = numChannels * sampleBitDepth / 8;
        if (bufferSize % frameSize != 0) {
            throwError("The buffer size must be a multiple of the frame size.",
                "Wave::readData(std::FILE*, uint32_t)");
        }

        int numFrames = numSamples / numChannels;
        uint32_t offset = 0;
        if (sampleBitDepth == 8) {
            while (offset < numSamples) {
                uint32_t numBytesToRead = std::min(bufferSize, dataChunkSize - offset);
                std::fread(pBuffer, 1, numBytesToRead, file);
                uint32_t numBufferFrames = numBytesToRead / frameSize;
                for (uint32_t y = 0; y < numBufferFrames; ++y) {
                    for (uint16_t x = 0; x < numChannels; ++x) {
                        uint32_t i = y * numChannels + x;
                        m_pData[offset + i] = static_cast<float>(pBuffer[i] / 127.5f - 1.f);
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
                        m_pData[offset + i] = static_cast<float>(pI16Buffer[i] / 32768.f);
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
                        m_pData[offset + i / 3] = static_cast<float>(
                            (pBuffer[i] << 8 | pBuffer[i + 1] << 16 | pBuffer[i + 2] << 24) / 2147483648.f);
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
                for (uint32_t y = 0; y < numBufferFrames; ++y) {
                    for (uint16_t x = 0; x < numChannels; ++x) {
                        uint32_t i = y * numChannels + x;
                        m_pData[offset + i] = static_cast<float>(pF32Buffer[i]);
                    }
                }
                offset += bufferSize / 4;
            }
        } else {
            std::free(pBuffer);
            throwError("Bit depth not supported.",
                "Wave::readData(std::FILE*, size_t)");
        }
        std::free(pBuffer);
    }

    /*!
     * \brief Fetches a buffer from the Wave object as a vector of floating-point values.
     * \param offset &mdash; the sample offset for the specified channel
     * \param sampleCount &mdash; the number of samples to fetch
     * \param channel &mdash; the index of the channel
     * \result An std::vector of floating-point values.
     */
    std::vector<float> Wave::getBuffer(uint32_t offset, uint32_t sampleCount, int channel) const
    {
        std::vector<float> buffer(sampleCount);
        uint16_t numChannels = int(m_waveProperties.getNumChannels());
        uint32_t absoluteOffset = numChannels * offset;
        for (uint32_t i = channel; i < sampleCount; i += numChannels) {
            buffer[i] = m_pData[absoluteOffset + i];
        }
        return buffer;
    }

    bool Wave::isEmpty() const
    {
        return (m_pData == nullptr);
    }

    float Wave::avgValue(int channel) const
    {
        uint16_t numChannels = m_waveProperties.getNumChannels();
        uint16_t offset = channel;
        double sum = 0.;
        for (uint32_t i = offset + numChannels; i < m_numSamples; i += numChannels) {
            sum += m_pData[i];
        }
        return float(sum / getNumFrames());
    }

    void Wave::changeVolume(float volume, int channel)
    { 
        uint16_t numChannels = m_waveProperties.getNumChannels();
        uint16_t offset = channel;
        if (getNumFrames() > 1) {
            float max = m_pData[offset];
            for (uint32_t i = offset + numChannels; i < m_numSamples; i += numChannels) {
                if (fabs(m_pData[i]) > max) {
                    max = fabs(m_pData[i]);
                }
            }
            float factor = volume / max;
            for (uint32_t i = offset; i < m_numSamples; i += numChannels) {
                m_pData[i] *= factor;
            }
        }
    }

    void Wave::downmixToMono()
    {
        uint32_t numSamples = getNumSamples();
        uint16_t numChannels = m_waveProperties.getNumChannels();
        uint32_t newDataChunkSize = m_waveProperties.getDataChunkSize() / numChannels;
        uint32_t newNumSamples = m_numSamples / numChannels;
        float* pNewData = reinterpret_cast<float*>(std::malloc(newNumSamples * sizeof(pNewData)));
        if (pNewData != nullptr) {
            for (uint32_t y = 0; y < numSamples; y += numChannels) {
                float sum = 0.f;
                for (uint16_t x = 0; x < numChannels; ++x) {
                    sum += m_pData[y + x];
                }
                pNewData[y / numChannels] = sum / float(numChannels);
            }
            std::free(m_pData);
            m_pData = pNewData;
            m_waveProperties.setBlockAlign(m_waveProperties.getBlockAlign() / numChannels);
            m_waveProperties.setNumBytesPerSecond(m_waveProperties.getNumBytesPerSecond() / numChannels);
            m_waveProperties.setRiffChunkSize(newDataChunkSize + 36);
            m_waveProperties.setDataChunkSize(newDataChunkSize);
            m_waveProperties.setNumChannels(1);
            m_numSamples = newNumSamples;
        }
    }

    Wave Wave::generateTestSineWave(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames)
    {
        Wave result(numFrames, 1);
        constexpr float coeff = 1.f - std::numeric_limits<float>::epsilon();
        float timeLimit = float(numFrames) / samplingFreq;
        float timeStep = timeLimit / (numFrames - 1);
        float omega = 2 * 3.1415927f * waveFreq;
        for (uint32_t i = 0; i < numFrames; ++i) {
            result.m_pData[i] = coeff * sinf(omega * (i * timeStep));
        }
        return result;
    }

    Wave Wave::generateTestSquareWave(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames)
    {
        Wave result(numFrames, 1);
        float timeStep = float(numFrames) / samplingFreq / (numFrames - 1);
        float omega = 2 * 3.1415927f * waveFreq;
        constexpr float coeff = 1.f - std::numeric_limits<float>::epsilon();
        for (uint32_t i = 0; i < numFrames; ++i) {
            result.m_pData[i] = sinf(omega * (i * timeStep)) < 0 ? -1.f * coeff : 1.f * coeff;
        }
        return result;
    }

    Wave Wave::generateTestTriangleWave(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames)
    {
        Wave result(numFrames, 1);
        float timeStep = float(numFrames) / samplingFreq / (numFrames - 1);
        constexpr float pi = 3.1415927f;
        constexpr float coeff = 2 * (1.f - std::numeric_limits<float>::epsilon()) / pi;
        float omega = 2 * pi * waveFreq;
        for (uint32_t i = 0; i < numFrames; ++i) {
            result.m_pData[i] = coeff * asinf(sinf(omega * (i * timeStep)));
        }
        return result;
    }

    /*!
     * \brief Checks if the Wave object is set to be stored using little-endian representation.
     * \result <b>true</b> if it is little-endian (RIFF), <b>false</b> if it is big-endian (RIFX).
     */
    bool Wave::isLittleEndian() const
    {
        return m_isLittleEndian;
    }

    float Wave::maxValue(int channel) const
    {
        uint16_t numChannels = m_waveProperties.getNumChannels();
        uint16_t offset = channel;
        float max = 0.f;
        if (getNumFrames() > 1) {
            max = m_pData[offset];
            for (uint32_t i = offset + numChannels; i < m_numSamples; i += numChannels) {
                if (m_pData[i] > max) {
                    max = m_pData[i];
                }
            }
        }
        return max;
    }

    float Wave::minValue(int channel) const
    {
        uint16_t numChannels = m_waveProperties.getNumChannels();
        uint16_t offset = channel;
        float min = 0.f;
        if (getNumFrames() > 1) {
            min = m_pData[offset];
            for (uint32_t i = offset + numChannels; i < m_numSamples; i += numChannels) {
                if (m_pData[i] < min) {
                    min = m_pData[i];
                }
            }
        }
        return min;
    }

    void Wave::reverse(int channel)
    {
        int numChannels = getNumChannels();
        int i = getNumSamples() - numChannels + channel;
        int j = channel;
        while (i > j) {
            float sample = m_pData[i];
            m_pData[i] = m_pData[j];
            m_pData[j] = sample;
            i -= numChannels;
            j += numChannels;
        }
    }

    void Wave::setLittleEndian(bool isLittleEndian)
    {
        m_isLittleEndian = isLittleEndian;
    }

    void Wave::swapChannels(int from, int to)
    {
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
            float sample = m_pData[i];
            m_pData[i] = m_pData[i + offset];
            m_pData[i + offset] = sample;
        }
    }

    void Wave::zeroInitHeader()
    {
        std::memset(&m_header, 0, sizeof(Header));
        std::memset(&m_dataSubChunk, 0, sizeof(DataSubChunk));
    }

    /*!
     * \brief Gets the number of bytes necessary to store the audio data with the current settings.
     * \result The data chunk size.
     */
    uint32_t Wave::getDataChunkSize() const
    {
        return m_waveProperties.getDataChunkSize();
    }

    /*!
     * \brief Gets the number of channels for the audio data.
     * \result The number of channels.
     */
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

    /*!
     * \brief Gets the number of bits currently used to represent a single audio sample.
     * \result The sample bit depth.
     */
    uint16_t Wave::getSampleBitDepth() const
    {
        return m_waveProperties.getNumBitsPerSample();
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

        bool isEndiannessMismatched = isCpuBigEndian() == m_isLittleEndian;
        uint32_t offset = 0;
        if (sampleBitDepth == 8) {
            while (offset < numSamples) {
                uint32_t numBytesToWrite = std::min(bufferSize, dataChunkSize - offset);
                uint32_t numBufferFrames = numBytesToWrite / frameSize;
                for (uint32_t y = 0; y < numBufferFrames; ++y) {
                    for (uint16_t x = 0; x < numChannels; ++x) {
                        int i = y * numChannels + x;
                        pBuffer[i] = uint8_t(roundf((m_pData[offset + i] + 1.f) * 127.5f));
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
                        pI16Buffer[i] = int16_t(m_pData[offset + i] * 32768.f);
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
                        int val = int(m_pData[offset + i / 3] * 2147483648.f);
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
                for (uint32_t y = 0; y < numBufferFrames; ++y) {
                    for (uint16_t x = 0; x < numChannels; ++x) {
                        uint32_t i = y * numChannels + x;
                        pF32Buffer[i] = m_pData[offset + i];
                    }
                }
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
        m_waveProperties.setNumBytesPerSecond(
            numChannels * m_waveProperties.getSamplingFrequency() * sampleBitDepth / 8);
        m_waveProperties.setBlockAlign(numChannels * sampleBitDepth / 8);
        m_waveProperties.setDataChunkSize(m_numSamples * sampleBitDepth / 8);
    }

    void Wave::setSampleRate(uint32_t sampleRate)
    {
        m_waveProperties.setSamplingFrequency(sampleRate);
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
