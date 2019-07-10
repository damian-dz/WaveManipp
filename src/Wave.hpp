#ifndef WAVE_H
#define WAVE_H

#include "WaveProperties.hpp"

namespace wm {
/*!
 * \brief A class that represents an audio file in memory.
 */
class Wave
{
    static_assert(std::numeric_limits<float>::is_iec559, "float must be IEEE 754");
    static_assert(sizeof(float) == sizeof(uint32_t), "float must be 32 bits wide");

private:
    struct Descriptor {
        char id[4];
        uint32_t size;
    };

    struct RIFFChunk {
        Descriptor descriptor;
        char type[4];
    };

    struct FmtSubChunk {
        Descriptor descriptor;
        uint16_t audioFormat;
        uint16_t numChannels;
        uint32_t sampleRate;
        uint32_t byteRate;
        uint16_t blockAlign;
        uint16_t bitsPerSample;
    };

    struct DataSubChunk {
        Descriptor descriptor;
    };

    struct Header {
        RIFFChunk riff;
        FmtSubChunk wave;
    };

    Header m_header;
    DataSubChunk m_dataSubChunk;
    bool m_isLittleEndian;
    uint32_t m_numSamples;
    float* m_pData;
    WaveProperties m_waveProperties;

    void copySamples(const float* source, float* destination, uint32_t count, uint32_t srcOffset = 0,
                     uint32_t destOffset = 0);
    void generateHeader();
    template <typename T> T reverseBytes(T val);
    void changeBufferEndianness(uint8_t* bytes, uint32_t sampleSize, uint32_t bufferSize);
    void findDataChunk(std::FILE* pFile);
    bool isCpuBigEndian() const;
    bool peekForId(const std::string& id, std::FILE* pFile);
    void setFourCharacterCodes();
    void setWaveProperties(bool onlyDataChunk = false);

public:
    Wave();
    Wave(const char* filename);
    Wave(const std::string& filename);
    Wave(uint32_t numFrames, uint16_t numChannels = 2, uint16_t bitDepth = 16, uint32_t sampleRate = 44100);
    Wave(const Wave& other);
    ~Wave();

    void open(const char* filename, uint32_t bufferSize = 24576);
    void open(const std::string& filename, uint32_t bufferSize = 24576);
    void readData(std::FILE* file, uint32_t bufferSize = 24576);

    Wave& append(const Wave& other);
    float avgValue(int channel = 0) const;
    void changeVolume(float volume, int channel = 0);
    void downmixToMono();
    static Wave generateSine(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames);
    static Wave generateSquare(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames);
    static Wave generateTriangle(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames);
    std::vector<float> getAveragedOutData(uint32_t binSize, int channel = 0) const;
    std::vector<float> getBuffer(uint32_t offset, uint32_t sampleCount, int channel = 0) const;
    bool isEmpty() const;
    bool isLittleEndian() const;
    bool isMono() const;
    float maxValue(int channel = 0) const;
    float minValue(int channel = 0) const;
    void reserveMemory(uint32_t numSamples, bool zeroInit = false);
    void resizeMemory(uint32_t numSamples, bool zeroInit = true);
    void reverse(int channel = 0);
    void setLittleEndian(bool isLittleEndian);
    void swapChannels(int from = 0, int to = 1);
    void zeroInitHeader();

    void saveAs(const char* filename, uint32_t bufferSize = 24576);
    void saveAs(const std::string& filename, uint32_t bufferSize = 24576);
    void writeData(std::FILE* file, uint32_t bufferSize = 24576);

    void setSampleBitDepth(uint16_t sampleBitDepth);
    void setSampleRate(uint32_t sampleRate);

    uint32_t getDataChunkSize() const;
    uint16_t getNumChannels() const;
    uint32_t getNumFrames() const;
    uint32_t getNumSamples() const;
    uint16_t getSampleBitDepth() const;

    friend std::ostream& operator<<(std::ostream& os, const Wave& wav);
};
}

#endif // !WAVE_H
