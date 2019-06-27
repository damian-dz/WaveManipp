#ifndef WAVE_H
#define WAVE_H

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

[[noreturn]] inline void throwError(const std::string &errorMsg, const std::string &funcName = "")
{
    std::string prefix = funcName != "" ? funcName + ": " : "";
    std::cerr << prefix << errorMsg << std::endl;
    throw errorMsg;
}

namespace wm {

    enum Channels
    {
        Left = 0, Right = 1, Both = 2,
    };

    /*!
     * \brief A class that represents an audio file in memory.
     */
    class Wave {
        static_assert(std::numeric_limits<float>::is_iec559, "float must be IEEE 754");
        static_assert(sizeof(float) == sizeof(uint32_t), "float must be 32 bits wide");

    private:
        struct chunk {
            char id[4];
            uint32_t size;
        };

        struct RIFFHeader {
            chunk descriptor;
            char type[4];
        };

        struct WAVEHeader {
            chunk descriptor;
            uint16_t audioFormat;
            uint16_t numChannels;
            uint32_t sampleRate;
            uint32_t byteRate;
            uint16_t blockAlign;
            uint16_t bitsPerSample;
        };

        struct DATAHeader {
            chunk descriptor;
        };

        struct CombinedHeader {
            RIFFHeader riff;
            WAVEHeader wave;
        };

        CombinedHeader m_combinedHeader;
        DATAHeader m_dataHeader;
        uint32_t m_dataSize;
        uint32_t m_numSamples;
        float *m_pData;

        void generateHeader(uint32_t dataLength, uint16_t numChannels, uint32_t sampleRate, uint16_t bitsPerSample);

    public:
        Wave();

        Wave(const char *c_pFilename);
        Wave(const std::string &filename);

        ~Wave();

        void open(const char *c_pFilename, size_t bufferSize = 24576);
        void open(const std::string &filename, size_t bufferSize = 24576);
        void readData(std::FILE *pFile, size_t bufferSize = 24576);

        std::vector<float> getBuffer(size_t offset, size_t sampleCount, Channels channels = Both) const;
        void changeVolume(float volume, Channels channels = Both);
        void downmixToMono();
        void swapChannels();
        void zeroInitHeader();

        void saveAs(const char *c_pFilename, uint32_t sampleRate = 44100, uint16_t sampleBitDepth = 16, size_t bufferSize = 24576);
        void saveAs(const std::string &filename, uint32_t sampleRate = 44100, uint16_t sampleBitDepth = 16, size_t bufferSize = 24576);
        void writeData(std::FILE *pFile, size_t bufferSize = 24576);

        uint32_t getDataSize() const;
        uint16_t getNumChannels() const;
        uint32_t getNumSamples() const;
        uint16_t getSampleBitDepth() const;
      
        friend std::ostream &operator<<(std::ostream &os, const Wave &wav);
    };
}

#endif // !WAVE_H
