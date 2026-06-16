#include "WaveMixer.hpp"


namespace wm {
WaveMixer::WaveMixer() :
    m_bitsPerSample(16),
    m_numChannels(1),
    m_sampleRate(default_sample_rate)
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

void WaveMixer::addTrack(uint16_t numChannels)
{
    Track track;
    track.setNumChannels(numChannels);
    m_tracks.push_back(Track());
}

void WaveMixer::addTrack(uint32_t offset, const Wave& wav)
{
    Track track;
    track.setNumChannels(wav.getNumChannels());
    track.addChunk(offset, wav);
    m_tracks.push_back(track);
}

void WaveMixer::addTrackAt(float timeInSec, const Wave& wav)
{
    uint32_t offset = static_cast<uint32_t>(std::round(static_cast<double>(timeInSec) * m_sampleRate));
    addTrack(offset, wav);
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

float WaveMixer::getTrackVolume(int trackIdx) const
{
    return m_tracks[trackIdx].getTrackVolume();
}

float WaveMixer::getTrackPan(int trackIdx) const
{
    return m_tracks[trackIdx].getTrackPan();
}

bool WaveMixer::isTrackMuted(int trackIdx) const
{
    return m_tracks[trackIdx].isTrackMuted();
}

bool WaveMixer::isTrackSolo(int trackIdx) const
{
    return m_tracks[trackIdx].isTrackSolo();
}

void WaveMixer::insertChunk(int trackIdx, uint32_t offset, const Wave& wav)
{
    m_tracks[trackIdx].addChunk(offset, wav);
}

void WaveMixer::removeTrack(int trackIdx)
{
    m_tracks.erase(m_tracks.begin() + trackIdx);
}

void WaveMixer::setBitsPerSample(uint16_t bitsPerSample)
{
    m_bitsPerSample = bitsPerSample;
}

void WaveMixer::setNumChannels(uint16_t numChannels)
{
    m_numChannels = numChannels;
}

void WaveMixer::setSampleRate(uint32_t sampleRate)
{
    m_sampleRate = sampleRate;
}

void WaveMixer::setTrackVolume(int trackIdx, float volume)
{
    m_tracks[trackIdx].setTrackVolume(volume);
}

void WaveMixer::setTrackPan(int trackIdx, float pan)
{
    m_tracks[trackIdx].setTrackPan(pan);
}

void WaveMixer::setTrackMuted(int trackIdx, bool muted)
{
    m_tracks[trackIdx].setTrackMuted(muted);
}

void WaveMixer::setTrackSolo(int trackIdx, bool solo)
{
    m_tracks[trackIdx].setTrackSolo(solo);
}

Wave WaveMixer::toWave_old()
{
    uint32_t numFrames = getNumFrames();
    Wave result(numFrames, m_numChannels, m_bitsPerSample, m_sampleRate);
    float* wave = result.audioData();
    std::memset(wave, 0, numFrames * uint64_t(m_numChannels) * sizeof(float));
    for (const Track& track : m_tracks) {
        uint32_t tStart = track.getMinStartOffset();
        uint32_t tEnd = track.getMaxEndOffset();

        uint32_t numTrackFrames = tEnd - tStart;
        uint32_t numTrackSamples = numTrackFrames * uint32_t(m_numChannels);
        float* trackData = reinterpret_cast<float*>(std::malloc(numTrackSamples * sizeof(float)));
        std::memset(trackData, 0, numTrackSamples * sizeof(float));
        uint16_t trackNumChannels = track.getNumChannels();
        float trackVolume = track.getTrackVolume();
        for (int i = 0; i < track.getNumChunks(); ++i) {
            Chunk chunk = track.getChunk(i);
            const float* chunkData = chunk.wave->constAudioData();
            uint32_t jStart = chunk.startOffset * trackNumChannels;
            uint32_t jEnd = chunk.endOffset * trackNumChannels;
            for (uint32_t j = jStart; j < jEnd; ++j) {
                trackData[j - tStart * trackNumChannels] += chunkData[j - jStart] * trackVolume;
            }
        }

        for (uint32_t y = tStart; y < tEnd; y++) {
            for (uint32_t x = 0; x < m_numChannels; x++) {
                int i = y * m_numChannels + x;
                int j = (y - tStart) * trackNumChannels + (x % trackNumChannels);
                wave[i] += trackData[j];
            }
        }
        std::free(trackData);
    }
    return result;
}

Wave WaveMixer::toWave()
{
    uint32_t numFrames = getNumFrames();
    Wave result(numFrames, m_numChannels, m_bitsPerSample, m_sampleRate);
    float* wave = result.audioData();
    std::memset(wave, 0, numFrames * uint64_t(m_numChannels) * sizeof(float));

    bool anySolo = false;
    for (const Track& track : m_tracks)
        if (track.isTrackSolo()) { anySolo = true; break; }

    for (const Track& track : m_tracks) {
        if (track.isTrackMuted()) continue;
        if (anySolo && !track.isTrackSolo()) continue;

        const float vol   = track.getTrackVolume();
        const float pan   = track.getTrackPan();
        const float gainL = vol * (pan <= 0.f ? 1.f : 1.f - pan);
        const float gainR = vol * (pan >= 0.f ? 1.f : 1.f + pan);
        const uint16_t srcCh = track.getNumChannels();

        for (int i = 0; i < track.getNumChunks(); ++i) {
            const Chunk chunk = track.getChunk(i);
            const float* src = chunk.wave->constAudioData();
            const uint32_t nf = chunk.endOffset - chunk.startOffset;
            for (uint32_t f = 0; f < nf; ++f) {
                const float sL = src[f * srcCh];
                const float sR = srcCh > 1 ? src[f * srcCh + 1] : sL;
                const uint32_t outF = chunk.startOffset + f;
                if (m_numChannels == 1) {
                    wave[outF] += (sL * gainL + sR * gainR) * 0.5f;
                } else {
                    wave[outF * m_numChannels + 0] += sL * gainL;
                    wave[outF * m_numChannels + 1] += sR * gainR;
                }
            }
        }
    }
    return result;
}

void WaveMixer::Track::addChunk(uint32_t startOffset, const Wave& wav)
{
    Chunk chunk;
    chunk.startOffset = startOffset;
    chunk.endOffset = startOffset + wav.getNumFrames();
    chunk.wave = &wav;
    m_chunks.push_back(chunk);
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

uint16_t WaveMixer::Track::getNumChannels() const
{
    return m_numChannels;
}

float WaveMixer::Track::getTrackVolume() const
{
    return m_trackVolume;
}

float WaveMixer::Track::getTrackPan() const
{
    return m_trackPan;
}

bool WaveMixer::Track::isTrackMuted() const
{
    return m_muted;
}

bool WaveMixer::Track::isTrackSolo() const
{
    return m_solo;
}

void WaveMixer::Track::setNumChannels(uint16_t numChannels)
{
    m_numChannels = numChannels;
}

void WaveMixer::Track::setTrackVolume(float volume)
{
    m_trackVolume = volume;
}

void WaveMixer::Track::setTrackPan(float pan)
{
    m_trackPan = pan;
}

void WaveMixer::Track::setTrackMuted(bool muted)
{
    m_muted = muted;
}

void WaveMixer::Track::setTrackSolo(bool solo)
{
    m_solo = solo;
}
}
