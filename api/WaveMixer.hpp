#pragma once

/*!
 * \file WaveMixer.hpp
 * \brief WaveMixer class for multi-track mixing of Wave chunks.
 */

#ifndef WAVE_MIXER_H
#define WAVE_MIXER_H

#include "Wave.hpp"

namespace wm {
/*!
 * \brief Mixes multiple Wave chunks across independent tracks.
 *
 * \details WaveMixer stores non-owning pointers to source Wave objects and renders
 * them into a new Wave when toWave() is called. Each track can contain one or more
 * chunks positioned by frame offset. Track volume is linear gain, pan ranges from
 * -1.0 (left) through 0.0 (center) to +1.0 (right), muted tracks are silent, and
 * soloed tracks suppress all non-soloed tracks.
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
        float getTrackPan() const;
        bool isTrackMuted() const;
        bool isTrackSolo() const;
        void setNumChannels(uint16_t numChannels);
        void setTrackVolume(float volume);
        void setTrackPan(float pan);
        void setTrackMuted(bool muted);
        void setTrackSolo(bool solo);

    private:
        std::vector<Chunk> m_chunks;
        uint16_t m_numChannels = 1;
        float m_trackVolume = 1.f;
        float m_trackPan    = 0.f;
        bool  m_muted       = false;
        bool  m_solo        = false;
    };

    uint16_t m_bitsPerSample;
    uint16_t m_numChannels;
    uint32_t m_sampleRate;
    std::vector<Track> m_tracks;

public:
    /*! \brief Creates an empty mono mixer at 16-bit, 44100 Hz output settings. */
    WAVEMANIPPAPI WaveMixer();
    /*! \brief Creates a mixer using a Wave's format and places it on the first track. */
    WAVEMANIPPAPI WaveMixer(const Wave& wav);
    /*! \brief Destroys the mixer without taking ownership of source waves. */
    WAVEMANIPPAPI ~WaveMixer();

    /*!
     * \brief Adds an empty track with the requested source channel count.
     * \param numChannels Number of channels in the track's source audio.
     */
    WAVEMANIPPAPI void addTrack(uint16_t numChannels);
    /*!
     * \brief Adds a track containing one Wave chunk at a frame offset.
     * \param offset Frame offset at which the chunk begins in the rendered output.
     * \param wav    Source wave to place on the new track.
     */
    WAVEMANIPPAPI void addTrack(uint32_t offset, const Wave& wav);
    /*!
     * \brief Adds a track containing one Wave chunk at a time offset in seconds.
     * \param timeInSec Time offset in seconds at which the chunk begins.
     * \param wav       Source wave to place on the new track.
     */
    WAVEMANIPPAPI void addTrackAt(float timeInSec, const Wave& wav);
    /*! \brief Returns the rendered length in frames. */
    WAVEMANIPPAPI uint32_t getNumFrames() const;
    /*! \brief Returns the number of tracks. */
    WAVEMANIPPAPI int getNumTracks() const;
    /*! \brief Returns a track's linear gain. \param trackIdx Zero-based track index. */
    WAVEMANIPPAPI float getTrackVolume(int trackIdx) const;
    /*! \brief Returns a track's pan value (-1 = full left, 0 = center, +1 = full right). \param trackIdx Zero-based track index. */
    WAVEMANIPPAPI float getTrackPan(int trackIdx) const;
    /*! \brief Returns true when a track is muted. \param trackIdx Zero-based track index. */
    WAVEMANIPPAPI bool  isTrackMuted(int trackIdx) const;
    /*! \brief Returns true when a track is soloed. \param trackIdx Zero-based track index. */
    WAVEMANIPPAPI bool  isTrackSolo(int trackIdx) const;
    /*!
     * \brief Inserts an additional Wave chunk into an existing track at a frame offset.
     * \param trackIdx Zero-based track index.
     * \param offset   Frame offset at which the chunk begins in the rendered output.
     * \param wav      Source wave to insert.
     */
    WAVEMANIPPAPI void insertChunk(int trackIdx, uint32_t offset, const Wave& wav);
    /*! \brief Removes a track and all its chunk references. \param trackIdx Zero-based track index. */
    WAVEMANIPPAPI void removeTrack(int trackIdx);
    /*! \brief Sets the output bit depth for rendered audio. */
    WAVEMANIPPAPI void setBitsPerSample(uint16_t bitsPerSample);
    /*! \brief Sets the output channel count. */
    WAVEMANIPPAPI void setNumChannels(uint16_t numChannels);
    /*! \brief Sets the output sample rate in Hz. */
    WAVEMANIPPAPI void setSampleRate(uint32_t sampleRate);
    /*!
     * \brief Sets a track's linear gain.
     * \param trackIdx Zero-based track index.
     * \param volume   Linear gain factor (1.0 = unity).
     */
    WAVEMANIPPAPI void setTrackVolume(int trackIdx, float volume);
    /*!
     * \brief Sets a track's pan value.
     * \param trackIdx Zero-based track index.
     * \param pan      Pan position: -1.0 = full left, 0.0 = center, +1.0 = full right.
     */
    WAVEMANIPPAPI void setTrackPan(int trackIdx, float pan);
    /*!
     * \brief Enables or disables track muting.
     * \param trackIdx Zero-based track index.
     * \param muted    Pass true to silence the track during rendering.
     */
    WAVEMANIPPAPI void setTrackMuted(int trackIdx, bool muted);
    /*!
     * \brief Enables or disables track soloing.
     * \param trackIdx Zero-based track index.
     * \param solo     Pass true to solo this track; non-soloed tracks are suppressed during rendering.
     */
    WAVEMANIPPAPI void setTrackSolo(int trackIdx, bool solo);
    /*! \brief Legacy renderer kept for compatibility. */
    WAVEMANIPPAPI Wave toWave_old();
    /*! \brief Renders all audible tracks into a new Wave. */
    WAVEMANIPPAPI Wave toWave();
};
}

#endif // !WAVE_MIXER_H
