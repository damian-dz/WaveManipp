#pragma once

#ifndef WAVE_MIXER_H
#define WAVE_MIXER_H

#include "Wave.hpp"

namespace wm {
/*!
 * \brief A class that makes it possibe to mix multiple audio tracks.
 */
class WaveMixer
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
    WAVEMANIPPAPI WaveMixer();
    WAVEMANIPPAPI WaveMixer(const Wave& wav);
    WAVEMANIPPAPI ~WaveMixer();

    WAVEMANIPPAPI void addTrack(uint16_t numChannels);
    WAVEMANIPPAPI void addTrack(uint32_t offset, const Wave& wav);
    WAVEMANIPPAPI void addTrackAt(float timeInSec, const Wave& wav);
    WAVEMANIPPAPI uint32_t getNumFrames() const;
    WAVEMANIPPAPI int getNumTracks() const;
    WAVEMANIPPAPI float getTrackVolume(int trackIdx) const;
    WAVEMANIPPAPI void insertChunk(int trackIdx, uint32_t offset, const Wave& wav);
    WAVEMANIPPAPI void removeTrack(int trackIdx);
    WAVEMANIPPAPI void setBitsPerSample(uint16_t bitsPerSample);
    WAVEMANIPPAPI void setNumChannels(uint16_t numChannels);
    WAVEMANIPPAPI void setSampleRate(uint32_t sampleRate);
    WAVEMANIPPAPI void setTrackVolume(int trackIdx, float volume);
    WAVEMANIPPAPI Wave toWave_old();
    WAVEMANIPPAPI Wave toWave();
};
}

#endif // !WAVE_MIXER_H
