#ifndef WAVEFILE_H
#define WAVEFILE_H

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <string>
#include <vector>

namespace wm {

    class WaveFile {
    private:
        struct chunk {
            char id[4];
            uint32_t size;
        };

        struct RIFFHeader {
            chunk descriptor;   // "RIFF"
            char type[4];       // "WAVE"
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

    public:
        WaveFile(uint8_t* data, uint32_t dataSize, uint16_t numChannels,
            uint32_t sampleRate, uint16_t bitsPerSample);
        WaveFile(uint32_t dataSize, uint16_t numChannels, uint32_t sampleRate,
            uint16_t bitsPerSample);
        WaveFile(const std::vector<float>& samples, uint16_t numChannels,
            uint32_t sampleRate, uint16_t bitsPerSample);
        WaveFile(const std::vector<float>& left, const std::vector<float>& right,
            uint32_t sampleRate, uint16_t bitsPerSample);
        WaveFile(const std::string fileName);
        WaveFile(const WaveFile& wav);

        ~WaveFile();

        void open(const std::string fileName);
        void saveAs(const std::string fileName);
        void changeVolume(float volume, int channelNumber = 0);
        std::vector<float> getAudioDataAsFloat(int channelNumber = 0);
        std::vector<int32_t> getAudioDataAsInt(int channelNumber = 0);
        void setAudioFromFloatData(const std::vector<float>& data, int channelNumber = 0);
        void setAudioFromIntData(const std::vector<int32_t>& data, int channelNumber = 0);
        std::vector<float> getAveragedOutAudioData(int binLength, int channelNumber = 0);
        void normalize(float newMin, float newMax, int channelNumber = 0);

        uint16_t numChannels() const;
        uint32_t sampleRate() const;
        void setSampleRate(int sampleRate);
        void setSampleBitDepth(int bitDepth);

        void append(const WaveFile& wav);
        void appendAudioFromFloatData(const std::vector<float>& data);
        void downmixToMono();

        static WaveFile mix(WaveFile& wav1, WaveFile& wav2, bool noClipping = true);
        static WaveFile join(WaveFile& wav1, WaveFile& wav2);

        friend std::ostream& operator<<(std::ostream& os, const WaveFile& wav);

    private:
        void generateHeader(uint32_t dataLength, uint16_t numChannels, uint32_t sampleRate, uint16_t bitsPerSample);

        CombinedHeader m_combinedHeader;
        DATAHeader m_dataHeader;
        uint32_t m_audioDataSize;
        uint8_t* m_pAudioData;
    };
}

#endif // WAVEFILE_H
