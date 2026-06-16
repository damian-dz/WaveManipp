#pragma once

/*!
 * \file Wave.hpp
 * \brief Core Wave class for loading, creating, and manipulating WAVE audio data.
 */

#ifndef WAVE_H
#define WAVE_H

#include "WaveProperties.hpp"
#include <memory>

namespace wm {
/*!
 * \brief Owns or shares decoded WAVE audio data in memory.
 *
 * \details Wave stores samples internally as normalized 32-bit floats while preserving
 * file-oriented metadata such as channel count, sample rate, bit depth, RIFF/RIFX
 * endianness, and chunk sizes. Samples are interleaved by frame: for stereo audio,
 * frame 0 is laid out as left sample 0, right sample 0, then frame 1, and so on.
 *
 * Copy construction and assignment are shallow and share the audio buffer. Mutating
 * methods detach from shared storage before writing.
 */
class WAVEMANIPPAPI Wave
{
    static_assert(std::numeric_limits<float>::is_iec559, "float must be IEEE 754");
    static_assert(sizeof(float) == sizeof(uint32_t), "float must be 32 bits wide");

private:
    struct Descriptor {
        char id[4];
        uint32_t size;
    };

    struct RIFFChunk {
        Descriptor descriptor;
        char type[4];
    };

    struct FmtSubChunk {
        Descriptor descriptor;
        uint16_t audioFormat;
        uint16_t numChannels;
        uint32_t sampleRate;
        uint32_t byteRate;
        uint16_t blockAlign;
        uint16_t bitsPerSample;
    };

    struct DataSubChunk {
        Descriptor descriptor;
    };

    struct Header {
        RIFFChunk riff;
        FmtSubChunk wave;
    };

    Header m_header;
    DataSubChunk m_dataSubChunk;
    bool m_isLittleEndian;
    uint32_t m_numSamples;
    std::shared_ptr<float[]> m_buffer;
    WaveProperties m_waveProperties;

    void detach();
    void copySamples(const float* source, float* destination, uint32_t count, uint32_t srcOffset = 0,
                     uint32_t destOffset = 0);
    void generateHeader();
    template <typename T> T reverseBytes(T val);
    void changeBufferEndianness(uint8_t* bytes, uint32_t sampleSize, uint32_t bufferSize);
    bool isCpuBigEndian() const;
    bool peekForId(const std::string& id, std::FILE* pFile);
    void setFourCharacterCodes();
    void setWaveProperties(bool onlyDataChunk = false);

public:
    /*! \brief Creates an empty wave with a zero-initialized header. */
    Wave();
    /*! \brief Loads a WAVE file into memory. \param filename Path to a RIFF or RIFX WAVE file. */
    Wave(const char* filename);
    /*! \brief Loads a WAVE file into memory. \param filename Path to a RIFF or RIFX WAVE file. */
    Wave(const std::string& filename);
    /*!
     * \brief Creates a silent wave with the requested format.
     * \param numFrames   Number of audio frames to allocate.
     * \param numChannels Number of interleaved channels.
     * \param bitDepth    Bit depth used when saving the file.
     * \param sampleRate  Sample rate in Hz.
     */
    Wave(uint32_t numFrames, uint16_t numChannels = 2, uint16_t bitDepth = 16, uint32_t sampleRate = default_sample_rate);
    /*! \brief Shares the audio buffer and metadata with another Wave. */
    Wave(const Wave& other);
    /*! \brief Releases this Wave's reference to the shared audio buffer. */
    ~Wave();

    /*!
     * \brief Opens and decodes a WAVE file.
     * \param filename   Path to a RIFF or RIFX WAVE file.
     * \param bufferSize Intermediate I/O buffer size in bytes; must align with the frame size.
     */
    void open(const char* filename, uint32_t bufferSize = default_buffer_size);
    /*!
     * \brief Opens and decodes a WAVE file.
     * \param filename   Path to a RIFF or RIFX WAVE file.
     * \param bufferSize Intermediate I/O buffer size in bytes; must align with the frame size.
     */
    void open(const std::string& filename, uint32_t bufferSize = default_buffer_size);
    /*!
     * \brief Reads audio data from an already-positioned file stream into this Wave.
     * \param file       Open file positioned at the data chunk payload.
     * \param bufferSize Intermediate read buffer size in bytes.
     */
    void readData(std::FILE* file, uint32_t bufferSize = default_buffer_size);

    /*! \brief Returns a read-only pointer to interleaved normalized float samples. */
    const float* constAudioData() const;
    /*! \brief Returns a mutable pointer to interleaved normalized float samples, detaching if shared. */
    float* audioData();

    /*! \brief Appends another Wave's audio after this Wave. \param other Source wave to append. */
    Wave& append(const Wave& other);
    /*!
     * \brief Returns the arithmetic mean sample value for one channel.
     * \param channel Zero-based channel index.
     */
    float avgValue(int channel = 0) const;
    /*!
     * \brief Multiplies every sample in one channel by a linear gain value.
     * \param volume  Linear gain factor (1.0 = unity, 0.5 = -6 dB, 2.0 = +6 dB).
     * \param channel Zero-based channel index.
     */
    void changeVolume(float volume, int channel = 0);
    /*! \brief Averages all channels into a single mono channel. */
    void downmixToMono();
    /*!
     * \brief Generates mono white noise in the range [-1, 1).
     * \param numFrames  Number of frames to generate.
     * \param bitDepth   Bit depth used when saving.
     * \param sampleRate Sample rate in Hz.
     */
    static Wave generateRandom(uint32_t numFrames, uint16_t bitDepth = 16, uint32_t sampleRate = default_sample_rate);
    /*!
     * \brief Generates mono white noise for a given duration.
     * \param duration   Duration in seconds.
     * \param bitDepth   Bit depth used when saving.
     * \param sampleRate Sample rate in Hz.
     */
    static Wave randomFromDuration(float duration, uint16_t bitDepth = 16, uint32_t sampleRate = default_sample_rate);
    /*!
     * \brief Generates a mono sine wave.
     * \param waveFreq      Frequency in Hz.
     * \param phaseShift    Phase offset in seconds.
     * \param numFrames     Number of frames to generate.
     * \param bitDepth      Bit depth used when saving.
     * \param sampleRate    Sample rate in Hz.
     * \param multiThreaded Enables OpenMP parallel generation when available.
     */
    static Wave generateSine(float waveFreq, float phaseShift,  uint32_t numFrames,
        uint16_t bitDepth = 16, uint32_t sampleRate = default_sample_rate, bool multiThreaded = false);
    /*!
     * \brief Generates a mono square wave.
     * \param waveFreq      Frequency in Hz.
     * \param phaseShift    Phase offset in seconds.
     * \param samplingFreq  Sample rate in Hz.
     * \param numFrames     Number of frames to generate.
     * \param multiThreaded Enables OpenMP parallel generation when available.
     */
    static Wave generateSquare(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames,
                               bool multiThreaded = false);
    /*!
     * \brief Generates a mono triangle wave.
     * \param waveFreq      Frequency in Hz.
     * \param phaseShift    Phase offset in seconds.
     * \param samplingFreq  Sample rate in Hz.
     * \param numFrames     Number of frames to generate.
     * \param multiThreaded Enables OpenMP parallel generation when available.
     */
    static Wave generateTriangle(float waveFreq, float phaseShift, uint32_t samplingFreq, uint32_t numFrames,
                                 bool multiThreaded = false);
    /*!
     * \brief Generates a mono click track at the specified tempo.
     * \param bpm         Tempo in beats per minute.
     * \param beatsPerBar Beats per bar (time signature numerator).
     * \param numFrames   Number of audio frames to generate.
     * \param bitDepth    Bit depth used when saving.
     * \param sampleRate  Sample rate in Hz.
     */
    static Wave generateClick(float bpm, uint8_t beatsPerBar, uint32_t numFrames,
                              uint16_t bitDepth = 16, uint32_t sampleRate = default_sample_rate);
    /*!
     * \brief Returns per-bin averaged channel data, optionally using absolute sample values.
     * \param binSize  Number of samples averaged per output bin.
     * \param absolute When true, uses absolute sample values before averaging.
     * \param channel  Zero-based channel index.
     */
    std::vector<float> getAveragedOutData(uint32_t binSize, bool absolute = false, int channel = 0) const;
    /*!
     * \brief Copies a contiguous channel segment into a vector.
     * \param offset      Starting sample offset within the channel.
     * \param sampleCount Number of samples to copy.
     * \param channel     Zero-based channel index.
     */
    std::vector<float> getBuffer(uint32_t offset, uint32_t sampleCount, int channel = 0) const;
    /*!
     * \brief Returns a downsampled view of a channel by averaging groups of samples.
     * \param offset               Starting sample offset within the channel.
     * \param squeezedSampleCount  Number of output samples.
     * \param squeezeFactor        Input-to-output decimation ratio.
     * \param absolute             When true, uses absolute sample values.
     * \param channel              Zero-based channel index.
     * \param multiThreaded        Enables OpenMP parallel processing when available.
     */
    std::vector<float> getSqueezedBuffer(uint32_t offset, uint32_t squeezedSampleCount, float squeezeFactor,
                                         bool absolute = false, int channel = 0, bool multiThreaded = false) const;
    /*!
     * \brief Returns an upsampled channel segment using simple linear stretching.
     * \param offset                Starting sample offset within the channel.
     * \param stretchedSampleCount  Number of output samples.
     * \param stretchFactor         Output-to-input expansion ratio.
     * \param channel               Zero-based channel index.
     */
    std::vector<float> getStretchedBuffer(uint32_t offset, uint32_t stretchedSampleCount, float stretchFactor,
                                          int channel = 0) const;
    /*!
     * \brief Returns the absolute peak sample value for one channel.
     * \param channel Zero-based channel index.
     */
    float getAbsPeak(int channel) const;
    /*!
     * \brief Inserts raw interleaved audio samples at a sample offset.
     * \param offset     Destination offset in interleaved samples (not frames).
     * \param audio      Pointer to interleaved float samples to copy.
     * \param numSamples Number of samples to insert.
     */
    void insertAudio(uint32_t offset, const float* audio, uint32_t numSamples);
    /*!
     * \brief Inserts raw interleaved audio samples at a sample offset.
     * \param offset Destination offset in interleaved samples (not frames).
     * \param audio  Vector of interleaved float samples to copy.
     */
    void insertAudio(uint32_t offset, std::vector<float>& audio);
    /*! \brief Returns true when this Wave has no allocated audio samples. */
    bool isEmpty() const;
    /*! \brief Returns true when this Wave will be serialized as RIFF little-endian data. */
    bool isLittleEndian() const;
    /*! \brief Returns true when the audio has one channel. */
    bool isMono() const;
    /*! \brief Returns the maximum sample value in one channel. \param channel Zero-based channel index. */
    float maxValue(int channel = 0) const;
    /*! \brief Returns the minimum sample value in one channel. \param channel Zero-based channel index. */
    float minValue(int channel = 0) const;
    /*!
     * \brief Allocates storage for a number of interleaved samples.
     * \param numSamples Total interleaved sample count to allocate.
     * \param zeroInit   When true, fills the new storage with zeros.
     */
    void reserveMemory(uint32_t numSamples, bool zeroInit = false);
    /*!
     * \brief Returns a linearly resampled copy at the requested sample rate.
     * \param targetSampleRate Output sample rate in Hz.
     */
    Wave resample(uint32_t targetSampleRate) const;
    /*!
     * \brief Resizes storage to a number of interleaved samples and updates chunk sizes.
     * \param numSamples Total interleaved sample count after resize.
     * \param zeroInit   When true, fills any newly added storage with zeros.
     */
    void resizeMemory(uint32_t numSamples, bool zeroInit = true);
    /*!
     * \brief Reverses the sample order of one channel in place.
     * \param channel Zero-based channel index.
     */
    void reverse(int channel = 0);
    /*!
     * \brief Replaces this Wave's interleaved sample data.
     * \param audio      Pointer to the new interleaved float samples.
     * \param numSamples Number of samples to copy.
     */
    void setAudio(const float* audio, uint32_t numSamples);
    /*!
     * \brief Replaces this Wave's interleaved sample data.
     * \param audio Vector of new interleaved float samples.
     */
    void setAudio(std::vector<float>& audio);
    /*!
     * \brief Selects RIFF little-endian or RIFX big-endian serialization.
     * \param isLittleEndian Pass true for RIFF (little-endian), false for RIFX (big-endian).
     */
    void setLittleEndian(bool isLittleEndian);
    /*!
     * \brief Swaps two channels in place.
     * \param from Zero-based index of the first channel.
     * \param to   Zero-based index of the second channel.
     */
    void swapChannels(int from = 0, int to = 1);
    /*! \brief Duplicates mono audio into a stereo buffer. */
    void upmixToStereo();
    /*! \brief Clears the internal RIFF and data chunk header structs. */
    void zeroInitHeader();

    /*!
     * \brief Saves this Wave as a RIFF/RIFX WAVE file.
     * \param filename   Destination file path.
     * \param bufferSize Intermediate I/O buffer size in bytes.
     */
    void saveAs(const char* filename, uint32_t bufferSize = default_buffer_size);
    /*!
     * \brief Saves this Wave as a RIFF/RIFX WAVE file.
     * \param filename   Destination file path.
     * \param bufferSize Intermediate I/O buffer size in bytes.
     */
    void saveAs(const std::string& filename, uint32_t bufferSize = default_buffer_size);
    /*!
     * \brief Writes the data chunk payload to an already-open file stream.
     * \param file       Open writable file positioned at the data chunk payload.
     * \param bufferSize Intermediate write buffer size in bytes.
     */
    void writeData(std::FILE* file, uint32_t bufferSize = default_buffer_size);

    /*!
     * \brief Sets the bit depth used when serializing samples.
     * \param sampleBitDepth Bit depth: 8, 16, 24, or 32.
     */
    void setSampleBitDepth(uint16_t sampleBitDepth);
    /*!
     * \brief Sets the sample rate metadata in Hz.
     * \param sampleRate Sample rate in Hz.
     */
    void setSampleRate(uint32_t sampleRate);

    /*! \brief Returns the byte size of the WAVE data chunk. */
    uint32_t getDataChunkSize() const;
    /*! \brief Returns the number of channels. */
    uint16_t getNumChannels() const;
    /*! \brief Returns the number of audio frames. */
    uint32_t getNumFrames() const;
    /*! \brief Returns the number of interleaved samples. */
    uint32_t getNumSamples() const;
    /*! \brief Returns the bit depth used for file serialization. */
    uint16_t getSampleBitDepth() const;
    /*! \brief Returns the sample rate in Hz. */
    uint32_t getSampleRate() const;

    /*! \brief Writes a human-readable summary of WAVE header metadata. */
    friend std::ostream& operator<<(std::ostream& os, const Wave& wav);

    /*!
     * \brief Accesses a mutable sample by frame index and channel index.
     * \param sample  Zero-based frame index.
     * \param channel Zero-based channel index.
     */
    float& operator()(uint32_t sample, int channel);
    /*!
     * \brief Accesses a read-only sample by frame index and channel index.
     * \param sample  Zero-based frame index.
     * \param channel Zero-based channel index.
     */
    const float& operator()(uint32_t sample, int channel) const;
    /*! \brief Shares another Wave's metadata and audio buffer with this Wave. */
    void operator=(const Wave& other);
};
}

#endif // !WAVE_H
