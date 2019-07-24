#include "WaveBuilder.hpp"

namespace wm {

WaveBuilder::WaveBuilder(const Wave& wav) :
    m_bitsPerSample(wav.getSampleBitDepth()),
    m_numChannels(wav.getNumChannels()),
    m_sampleRate(wav.getSampleRate())
{
    m_wavPtrs.push_back(&wav);
}

WaveBuilder::WaveBuilder(const WaveBuilder& other) :
    m_wavPtrs(other.m_wavPtrs),
    m_bitsPerSample(other.m_bitsPerSample),
    m_numChannels(other.m_numChannels),
    m_sampleRate(other.m_sampleRate)
{

}

WaveBuilder::~WaveBuilder()
{

} 

WaveBuilder& WaveBuilder::append(const Wave& wav)
{
    if (wav.getSampleBitDepth() != m_bitsPerSample ||
        wav.getSampleRate() != m_sampleRate ||
        wav.getNumChannels() != m_numChannels) {
        throwError("Wave parameters mismatched.",
                    "WaveBuilder::append(const Wave&)");
    }
    m_wavPtrs.push_back(&wav);
    return *this;
}

void WaveBuilder::insert(int idx, const Wave& wav)
{
    if (wav.getSampleBitDepth() != m_bitsPerSample ||
        wav.getSampleRate() != m_sampleRate ||
        wav.getNumChannels() != m_numChannels) {
        throwError("Wave parameters mismatched.",
                   "WaveBuilder::insert(const Wave&)");
    }
    m_wavPtrs.insert(m_wavPtrs.begin() + idx, &wav);
}

Wave WaveBuilder::toWave()
{
    uint32_t totalNumSamples = 0;
    for (size_t i = 0; i < m_wavPtrs.size(); ++i) {
        totalNumSamples += m_wavPtrs[i]->getNumSamples();
    }
    uint32_t offset = 0;
    Wave result(totalNumSamples / m_numChannels, m_numChannels, m_bitsPerSample, m_sampleRate);
    for (size_t i = 0; i < m_wavPtrs.size(); ++i) {
        uint32_t numSamples = m_wavPtrs[i]->getNumSamples();
        result.insertAudio(offset, m_wavPtrs[i]->constAudioData(), m_wavPtrs[i]->getNumSamples());
        offset += numSamples;
    }
    return result;
}
}