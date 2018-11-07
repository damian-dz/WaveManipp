#include "WaveBuilder.h"

#include <iostream>

namespace wm {


    WaveBuilder::WaveBuilder(const WaveFile &wav) :
        m_pFinalData(NULL)
    {
        m_bitsPerSample = wav.sampleBitDepth();
        m_numChannels = wav.numChannels();
        m_sampleRate = wav.sampleRate();
        m_pointers.push_back(wav.audioData());
        m_sizes.push_back(wav.audioDataSize());
    }

    WaveBuilder::~WaveBuilder()
    {
        std::cout << "WaveBuilder destroyed" << std::endl;
     //   clearAudioData();
    }

    void WaveBuilder::clearAudioData()
    {
        if (m_pFinalData != NULL) {
            std::free(m_pFinalData);
            m_pFinalData = NULL;
            std::cout << "WaveBuilder data freed" << std::endl;
        }
    }

    void WaveBuilder::append(const WaveFile &wav)
    {
        if (wav.sampleBitDepth() != m_bitsPerSample ||
            wav.sampleRate() != m_sampleRate ||
            wav.numChannels() != m_numChannels) {
            throw "WaveFile parameters mismatched";
        }
        m_pointers.push_back(wav.audioData());
        m_sizes.push_back(wav.audioDataSize());
    }

    WaveFile WaveBuilder::toWaveFile()
    {
        uint32_t dataSize = 0;
        for (uint32_t ds : m_sizes) {
            dataSize += ds;
        }
        m_pFinalData = (uint8_t*)std::malloc(dataSize);
        for (int i = 0; i < m_pointers.size(); ++i) {
            std::memcpy(&m_pFinalData[0], m_pointers[i], m_sizes[i]);
            m_pFinalData += m_sizes[i];
        }
        m_pFinalData -= dataSize;
        return WaveFile(m_pFinalData, dataSize, m_numChannels, m_sampleRate, m_bitsPerSample);
    }

}