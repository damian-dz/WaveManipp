#pragma once

#include "WaveFile.h"

#include <vector>

namespace wm {

    class WaveBuilder {
    public:
        WaveBuilder(const WaveFile &wav);
        ~WaveBuilder();

        void clearAudioData();
        void append(const WaveFile &wav);
        WaveFile toWaveFile();

    private:
        std::vector<uint8_t*> m_pointers;
        std::vector<uint32_t> m_sizes;
        uint16_t m_numChannels;
        uint32_t m_sampleRate;
        uint16_t m_bitsPerSample;

        uint8_t *m_pFinalData;
    };

}
