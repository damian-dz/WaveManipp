#pragma once

#ifndef WAVE_MIXER_H
#define WAVE_MIXER_H

#include "Wave.hpp"

namespace wm {
/*!
 * \brief A class that makes it possibe to mix multiple audio tracks.
 */
class WAVEMANIPPAPI WaveMixer
{
    struct Chunk
    {
        uint32_t startOffset = 0;
        uint32_t endOffset = 0;
        const Wave* wave = nullptr;
    };

    struct Track
    {
        void addChunk(uint32_t startOffset, const Wave& wav);
        Chunk getChunk(int idx) const;
        Chunk getLastChunk() const;
        uint32_t getMaxEndOffset() const;
        uint32_t getMinStartOffset() const;
        int getNumChunks() const;
        uint16_t getNumChannels() const;
        float getTrackVolume() const;
        void setNumChannels(uint16_t numChannels);
        void setTrackVolume(float volume);

    private:
        std::vector<Chunk> m_chunks;
        uint16_t m_numChannels;
        float m_trackVolume = 1.f;
    };

    uint16_t m_bitsPerSample;
    uint16_t m_numChannels;
    uint32_t m_sampleRate;
    std::vector<Track> m_tracks;

public:
    WaveMixer();
    WaveMixer(const Wave& wav);
    ~WaveMixer();

    void addTrack(uint16_t numChannels);
    void addTrack(uint32_t offset, const Wave& wav);
    void addTrackAt(float timeInSec, const Wave& wav);
    uint32_t getNumFrames() const;
    int getNumTracks() const;
    float getTrackVolume(int trackIdx) const;
    void insertChunk(int trackIdx, uint32_t offset, const Wave& wav);
    void removeTrack(int trackIdx);
    void setBitsPerSample(uint16_t bitsPerSample);
    void setNumChannels(uint16_t numChannels);
    void setSampleRate(uint32_t sampleRate);
    void setTrackVolume(int trackIdx, float volume);
    Wave toWave_old();
    Wave toWave();
};
}

#endif // !WAVE_MIXER_H
