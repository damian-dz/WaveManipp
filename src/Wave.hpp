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
        uint32_t m_numSamples;
        float* m_pData;
        bool m_isLittleEndian;

        WaveProperties m_waveProperties;

        void generateHeader();
        template <typename T> T reverseBytes(T val);
        void changeBufferEndianness(uint8_t* bytes, size_t sampleLength, size_t bufferLength);
        bool isCpuBigEndian() const;
        void setWaveProperties(bool onlyDataChunk = false);

        void findDataChunk(std::FILE* pFile);
        bool peekForId(const std::string& id, std::FILE* pFile);


    public:
        Wave();

        Wave(const char* filename);
        Wave(const std::string& filename);

        ~Wave();

        void open(const char* filename, uint32_t bufferSize = 24576);
        void open(const std::string& filename, uint32_t bufferSize = 24576);
        void readData(std::FILE* file, uint32_t bufferSize = 24576);

        std::vector<float> getBuffer(uint32_t offset, uint32_t sampleCount, int channel = 0) const;
        void changeVolume(float volume, int channel = 0);
        void downmixToMono();
        bool isLittleEndian() const;
        void reverse(int channel = 0);
        void setLittleEndian(bool val);
        void swapChannels();
        void zeroInitHeader();

        void saveAs(const char* filename, uint32_t bufferSize = 24576);
        void saveAs(const std::string& filename, uint32_t bufferSize = 24576);
        void writeData(std::FILE* file, uint32_t bufferSize = 24576);

        void setSampleBitDepth(uint16_t sampleBitDepth);
        void setSampleRate(uint32_t sampleRate);

        uint32_t getDataChunkSize() const;
        uint16_t getNumChannels() const;
        uint32_t getNumSamples() const;
        uint16_t getSampleBitDepth() const;

        friend std::ostream& operator<<(std::ostream& os, const Wave& wav);
    };

}

#endif // !WAVE_H
