#include "WaveMixer.hpp"


namespace wm {

WaveMixer::WaveMixer() :
    m_bitsPerSample(16),
    m_numChannels(1),
    m_sampleRate(44100)
{
}

WaveMixer::WaveMixer(const Wave& wav) :
    m_bitsPerSample(wav.getSampleBitDepth()),
    m_numChannels(wav.getNumChannels()),
    m_sampleRate(wav.getSampleRate())
{
    Track track;
    track.addChunk(0, wav);
    m_tracks.push_back(track);
}

WaveMixer::~WaveMixer()
{
}

void WaveMixer::addTrack()
{
    m_tracks.push_back(Track());
}

uint32_t WaveMixer::getNumFrames() const
{
    uint32_t max = 0;
    for (const Track &track : m_tracks) {
        uint32_t endOffset = track.getMaxEndOffset();
        if (endOffset > max) {
            max = endOffset;
        }
    }
    return max;
}

int WaveMixer::getNumTracks() const
{
    return int(m_tracks.size());
}

void WaveMixer::insertWave(int trackIdx, uint32_t offset, const Wave& wav)
{
    m_tracks[trackIdx].addChunk(offset, wav);
}

void WaveMixer::setTrackVolume(int trackIdx, float volume)
{
    m_tracks[trackIdx].setTrackVolume(volume);
}

Wave WaveMixer::toWave()
{
    uint32_t numFrames = getNumFrames();
    Wave result(numFrames, m_numChannels, m_bitsPerSample, m_sampleRate);
    float* data = result.audioData();
    std::memset(data, 0, numFrames * uint64_t(m_numChannels) * sizeof(float));
    for (const Track& track : m_tracks) {
        uint32_t tStart = track.getMinStartOffset();
        uint32_t tEnd = track.getMaxEndOffset();
        uint32_t numTrackFrames = tEnd - tStart;
        uint32_t numTrackSamples = numTrackFrames * uint32_t(m_numChannels);
        float* trackData = reinterpret_cast<float*>(std::malloc(numTrackSamples * sizeof(float)));
        std::memset(trackData, 0, numTrackSamples * sizeof(float));
        float trackVolume = track.getTrackVolume();
        for (int i = 0; i < track.getNumChunks(); ++i) {
            Chunk chunk = track.getChunk(i);
            const float* chunkData = chunk.data->constAudioData();
            uint32_t jStart = chunk.startOffset * m_numChannels;
            uint32_t jEnd = chunk.endOffset * m_numChannels;
            for (int j = jStart; j < jEnd; ++j) {
                trackData[j - jStart] = chunkData[j - jStart] * trackVolume;
            }
        }
        for (int i = 0; i < numTrackSamples; ++i) {
            data[i + tStart] += trackData[i];
        }
        std::free(trackData);
    }
    return result;
}

void WaveMixer::Track::addChunk(uint32_t startOffset, uint32_t endOffset, const Wave* wavData)
{
    Chunk chunk;
    chunk.startOffset = startOffset;
    chunk.endOffset = endOffset;
    chunk.data = wavData;
    m_chunks.push_back(chunk);
}

void WaveMixer::Track::addChunk(uint32_t startOffset, const Wave& wav)
{
    addChunk(startOffset, startOffset + wav.getNumFrames(), &wav);
}

WaveMixer::Chunk WaveMixer::Track::getChunk(int idx) const
{
    return m_chunks[idx];
}

WaveMixer::Chunk WaveMixer::Track::getLastChunk() const
{
    return m_chunks[m_chunks.size() - 1];
}

uint32_t WaveMixer::Track::getMaxEndOffset() const
{
    uint32_t maxEndOffset = 0;
    for (const Chunk& chunk : m_chunks) {
        if (chunk.endOffset > maxEndOffset) {
            maxEndOffset = chunk.endOffset;
        }
    }
    return maxEndOffset;
}

uint32_t WaveMixer::Track::getMinStartOffset() const
{
    if (m_chunks.size() < 1) {
        return 0;
    }
    uint32_t minStartOffset = m_chunks[0].startOffset;
    for (size_t i = 1; i < m_chunks.size(); ++i) {
        if (m_chunks[i].startOffset < minStartOffset) {
            minStartOffset = m_chunks[i].startOffset;
        }
    }
    return minStartOffset;
}

int WaveMixer::Track::getNumChunks() const
{
    return int(m_chunks.size());
}

float WaveMixer::Track::getTrackVolume() const
{
    return m_trackVolume;
}

void WaveMixer::Track::setTrackVolume(float volume)
{
    m_trackVolume = volume;
}
}
