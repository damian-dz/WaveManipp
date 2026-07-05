#include <doctest/doctest.h>

#include "PcmCodec.hpp"
#include "WaveFileReader.hpp"
#include "WaveFileFormat.hpp"
#include "Wave.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

// docs/AUDIO_STREAMING_PLAN.md Milestone 1 item 4: wm::WaveFileReader parses WAVE
// metadata without decoding the whole data chunk, then converts bounded frame ranges
// to normalized interleaved floats on demand. Wave and WaveFileReader share the private
// header parser and integer PCM codec, so handcrafted fixtures and direct codec tests
// assert independently-known metadata, integer codes, byte order, and sample values;
// Wave comparisons remain useful integration-consistency checks.

namespace {

std::filesystem::path tempWavPath(const char* suffix)
{
    return std::filesystem::temp_directory_path() / (std::string("wavemanipp_wfr_test_") + suffix + ".wav");
}

void writeU16(std::ofstream& out, uint16_t value, bool littleEndian)
{
    const uint8_t bytes[2] = {
        uint8_t(littleEndian ? value : value >> 8),
        uint8_t(littleEndian ? value >> 8 : value)
    };
    out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void writeU32(std::ofstream& out, uint32_t value, bool littleEndian)
{
    const uint8_t bytes[4] = {
        uint8_t(littleEndian ? value : value >> 24),
        uint8_t(littleEndian ? value >> 8 : value >> 16),
        uint8_t(littleEndian ? value >> 16 : value >> 8),
        uint8_t(littleEndian ? value >> 24 : value)
    };
    out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void writePcm24Fixture(const std::filesystem::path& path, bool littleEndian,
                       const std::vector<int32_t>& samples)
{
    const uint32_t dataSize = uint32_t(samples.size() * 3);
    std::ofstream out(path, std::ios::binary);
    out.write(littleEndian ? "RIFF" : "RIFX", 4);
    writeU32(out, 36 + dataSize, littleEndian);
    out.write("WAVEfmt ", 8);
    writeU32(out, 16, littleEndian);
    writeU16(out, 1, littleEndian);
    writeU16(out, 1, littleEndian);
    writeU32(out, 48000, littleEndian);
    writeU32(out, 48000 * 3, littleEndian);
    writeU16(out, 3, littleEndian);
    writeU16(out, 24, littleEndian);
    out.write("data", 4);
    writeU32(out, dataSize, littleEndian);

    for (int32_t sample : samples) {
        const uint32_t bits = static_cast<uint32_t>(sample) & 0x00FFFFFFu;
        const uint8_t bytes[3] = {
            uint8_t(littleEndian ? bits : bits >> 16),
            uint8_t(bits >> 8),
            uint8_t(littleEndian ? bits >> 16 : bits)
        };
        out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
    }
}

void writeExtensibleFixture(const std::filesystem::path& path, uint16_t subFormat,
                            const std::vector<uint8_t>& sampleBytes, uint16_t bitsPerSample)
{
    constexpr bool littleEndian = true;
    constexpr uint16_t channels = 1;
    constexpr uint32_t sampleRate = 48000;
    constexpr uint32_t fmtChunkSize = 40;
    const uint16_t blockAlign = channels * bitsPerSample / 8;
    const uint32_t dataSize = uint32_t(sampleBytes.size());

    std::ofstream out(path, std::ios::binary);
    out.write("RIFF", 4);
    writeU32(out, 4 + 8 + fmtChunkSize + 8 + dataSize, littleEndian);
    out.write("WAVEfmt ", 8);
    writeU32(out, fmtChunkSize, littleEndian);
    writeU16(out, 0xFFFE, littleEndian);
    writeU16(out, channels, littleEndian);
    writeU32(out, sampleRate, littleEndian);
    writeU32(out, sampleRate * blockAlign, littleEndian);
    writeU16(out, blockAlign, littleEndian);
    writeU16(out, bitsPerSample, littleEndian);
    writeU16(out, 22, littleEndian); // extension size
    writeU16(out, bitsPerSample, littleEndian); // valid bits per sample
    writeU32(out, 0, littleEndian); // channel mask
    writeU32(out, subFormat, littleEndian); // PCM or IEEE-float GUID data1
    const uint8_t guidTail[12] = {0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xAA, 0x00, 0x38, 0x9B, 0x71};
    out.write(reinterpret_cast<const char*>(guidTail), sizeof(guidTail));
    out.write("data", 4);
    writeU32(out, dataSize, littleEndian);
    out.write(reinterpret_cast<const char*>(sampleBytes.data()), sampleBytes.size());
}

wm::Wave makeTestWave(uint32_t numFrames, uint16_t channels, uint16_t bitDepth, uint32_t sampleRate,
                      bool littleEndian)
{
    wm::Wave wave(numFrames, channels, bitDepth, sampleRate);
    std::vector<float> samples(size_t(numFrames) * channels);
    for (size_t i = 0; i < samples.size(); ++i) {
        // Deterministic, non-trivial pattern that sweeps close to both extremes
        // (quantization edge cases for 8/16/24-bit) without relying on randomness.
        const float t = float(i) / float(samples.size());
        samples[i] = std::sin(37.0f * t) * (0.999f - 0.001f * float(i % 7));
    }
    wave.setAudio(samples);
    wave.setLittleEndian(littleEndian);
    return wave;
}

void checkMetadataMatches(const wm::Wave& oracle, const wm::WaveFileReader& reader)
{
    CHECK(reader.getNumChannels() == oracle.getNumChannels());
    CHECK(reader.getSampleRate() == oracle.getSampleRate());
    CHECK(reader.getSampleBitDepth() == oracle.getSampleBitDepth());
    CHECK(reader.getNumFrames() == oracle.getNumFrames());
    CHECK(reader.isLittleEndian() == oracle.isLittleEndian());
}

void checkSequentialReadMatchesOracle(const std::filesystem::path& path)
{
    wm::Wave oracle(path.string());
    wm::WaveFileReader reader(path.string());
    checkMetadataMatches(oracle, reader);

    const uint16_t channels = oracle.getNumChannels();
    const uint32_t numFrames = oracle.getNumFrames();
    std::vector<float> buffer(size_t(numFrames) * channels, -2.f);

    const uint32_t got = reader.readFrames(0, numFrames, buffer.data());
    REQUIRE(got == numFrames);

    const float* expected = oracle.constAudioData();
    for (size_t i = 0; i < buffer.size(); ++i) {
        CHECK(buffer[i] == expected[i]);
    }
}

template <typename Sample, unsigned Bits>
void checkPcmQuantizer()
{
    constexpr int64_t scale = int64_t{1} << (Bits - 1);
    constexpr int32_t minimum = static_cast<int32_t>(-scale);
    constexpr int32_t maximum = static_cast<int32_t>(scale - 1);
    const Sample lsb = Sample{1} / static_cast<Sample>(scale);
    const auto quantize = [](Sample value) { return wm::detail::quantizeSignedPcm<Sample, Bits>(value); };
    const auto dequantize = [](int32_t value) { return wm::detail::dequantizeSignedPcm<Sample, Bits>(value); };

    CAPTURE(sizeof(Sample));
    CAPTURE(Bits);

    CHECK(quantize(Sample{0}) == 0);
    CHECK(quantize(static_cast<Sample>(0.4) * lsb) == 0);
    CHECK(quantize(static_cast<Sample>(-0.4) * lsb) == 0);
    CHECK(quantize(static_cast<Sample>(0.5) * lsb) == 1);
    CHECK(quantize(static_cast<Sample>(-0.5) * lsb) == -1);
    CHECK(quantize(static_cast<Sample>(0.6) * lsb) == 1);
    CHECK(quantize(static_cast<Sample>(-0.6) * lsb) == -1);
    CHECK(quantize(Sample{-1}) == minimum);
    CHECK(quantize(Sample{1}) == maximum);
    CHECK(quantize(Sample{-2}) == minimum);
    CHECK(quantize(Sample{2}) == maximum);
    CHECK(quantize(std::numeric_limits<Sample>::quiet_NaN()) == 0);
    CHECK(quantize(-std::numeric_limits<Sample>::infinity()) == minimum);
    CHECK(quantize(std::numeric_limits<Sample>::infinity()) == maximum);

    CHECK(dequantize(minimum) == Sample{-1});
    CHECK(dequantize(0) == Sample{0});

    if constexpr (Bits <= std::numeric_limits<Sample>::digits) {
        for (int32_t code : {minimum, int32_t(minimum + 1), int32_t(-1), int32_t(0), int32_t(1),
                             int32_t(maximum - 1), maximum}) {
            const Sample decoded = dequantize(code);
            CHECK(quantize(decoded) == code);
        }
    }
}

} // namespace

// ============================================================================
// Format/endianness sweep: sequential full-range read vs wm::Wave::open()
// ============================================================================

TEST_CASE("WaveFileReader matches wm::Wave::open() across every supported bit depth and endianness")
{
    for (uint16_t bitDepth : {8, 16, 24, 32}) {
        for (bool littleEndian : {true, false}) {
            for (uint16_t channels : {1, 2}) {
                CAPTURE(bitDepth);
                CAPTURE(littleEndian);
                CAPTURE(channels);

                wm::Wave wave = makeTestWave(500, channels, bitDepth, 44100, littleEndian);
                const auto path = tempWavPath("sweep");
                wave.saveAs(path.string());

                checkSequentialReadMatchesOracle(path);

                std::filesystem::remove(path);
            }
        }
    }
}

TEST_CASE("internal file positioning preserves offsets beyond 2 GiB")
{
    wm::detail::FileHandle file(std::tmpfile());
    REQUIRE(static_cast<bool>(file));

    constexpr uint64_t offset = (uint64_t(1) << 31) + 123;
    REQUIRE(wm::detail::seekFile(file.get(), offset));

    uint64_t actualOffset = 0;
    REQUIRE(wm::detail::tellFile(file.get(), actualOffset));
    CHECK(actualOffset == offset);
}

TEST_CASE("integer PCM quantization is consistent for float and double samples")
{
    checkPcmQuantizer<float, 8>();
    checkPcmQuantizer<float, 16>();
    checkPcmQuantizer<float, 24>();
    checkPcmQuantizer<float, 32>();
    checkPcmQuantizer<double, 8>();
    checkPcmQuantizer<double, 16>();
    checkPcmQuantizer<double, 24>();
    checkPcmQuantizer<double, 32>();
}

TEST_CASE("unsigned 8-bit WAVE PCM uses 128 as digital zero")
{
    CHECK(wm::detail::encodePcm8Sample(-1.f) == 0);
    CHECK(wm::detail::encodePcm8Sample(0.f) == 128);
    CHECK(wm::detail::encodePcm8Sample(1.f) == 255);
    CHECK(wm::detail::decodePcm8Sample<float>(0) == -1.f);
    CHECK(wm::detail::decodePcm8Sample<float>(128) == 0.f);
    CHECK(wm::detail::decodePcm8Sample<float>(255) == 127.f / 128.f);

    for (uint32_t byte = 0; byte <= 255; ++byte) {
        const float decoded = wm::detail::decodePcm8Sample<float>(static_cast<uint8_t>(byte));
        CHECK(wm::detail::encodePcm8Sample(decoded) == byte);
    }
}

TEST_CASE("24-bit PCM decoding handles signed extrema in RIFF and RIFX files")
{
    const std::vector<int32_t> samples = {-8388608, -8388607, -1, 0, 1, 8388607};

    for (bool littleEndian : {true, false}) {
        CAPTURE(littleEndian);
        const auto path = tempWavPath(littleEndian ? "pcm24_riff" : "pcm24_rifx");
        writePcm24Fixture(path, littleEndian, samples);

        wm::WaveFileReader reader(path.string());
        REQUIRE(reader.getNumFrames() == samples.size());
        CHECK(reader.isLittleEndian() == littleEndian);

        std::vector<float> decoded(samples.size());
        REQUIRE(reader.readFrames(0, uint32_t(samples.size()), decoded.data()) == samples.size());
        for (size_t i = 0; i < samples.size(); ++i) {
            CHECK(decoded[i] == float(samples[i]) / 8388608.f);
        }

        reader.close();
        std::filesystem::remove(path);
    }
}

TEST_CASE("24-bit PCM encoding rounds to nearest and writes file byte order")
{
    struct EncodeCase {
        float input;
        int32_t expected;
    };

    constexpr float lsb = 1.f / 8388608.f;
    const std::vector<EncodeCase> cases = {
        {0.f, 0},
        {0.4f * lsb, 0},
        {0.5f * lsb, 1},
        {0.6f * lsb, 1},
        {-0.4f * lsb, 0},
        {-0.5f * lsb, -1},
        {-0.6f * lsb, -1},
        {-1.f, -8388608},
        {1.f, 8388607},
        {-2.f, -8388608},
        {2.f, 8388607},
        {std::numeric_limits<float>::quiet_NaN(), 0}
    };

    for (bool littleEndian : {true, false}) {
        CAPTURE(littleEndian);
        for (const EncodeCase& testCase : cases) {
            CAPTURE(testCase.input);
            uint8_t encoded[3]{};
            wm::detail::encodePcm24Sample(testCase.input, littleEndian, encoded);

            const uint32_t bits = static_cast<uint32_t>(testCase.expected) & 0x00FFFFFFu;
            CHECK(encoded[0] == uint8_t(littleEndian ? bits : bits >> 16));
            CHECK(encoded[1] == uint8_t(bits >> 8));
            CHECK(encoded[2] == uint8_t(littleEndian ? bits >> 16 : bits));
        }
    }
}

TEST_CASE("WAVE_FORMAT_EXTENSIBLE resolves PCM and IEEE-float sub-formats")
{
    SUBCASE("16-bit PCM")
    {
        const auto path = tempWavPath("extensible_pcm");
        const std::vector<int16_t> samples = {-32768, -1, 0, 1, 32767};
        std::vector<uint8_t> bytes;
        bytes.reserve(samples.size() * 2);
        for (int16_t sample : samples) {
            const uint16_t bits = static_cast<uint16_t>(sample);
            bytes.push_back(uint8_t(bits));
            bytes.push_back(uint8_t(bits >> 8));
        }
        writeExtensibleFixture(path, 1, bytes, 16);

        wm::WaveFileReader reader(path.string());
        CHECK(reader.getAudioFormat() == 1);
        CHECK(reader.getSampleBitDepth() == 16);
        REQUIRE(reader.getNumFrames() == samples.size());

        std::vector<float> decoded(samples.size());
        REQUIRE(reader.readFrames(0, uint32_t(samples.size()), decoded.data()) == samples.size());
        for (size_t i = 0; i < samples.size(); ++i) {
            CHECK(decoded[i] == float(samples[i]) / 32768.f);
        }

        reader.close();
        std::filesystem::remove(path);
    }

    SUBCASE("32-bit IEEE float")
    {
        const auto path = tempWavPath("extensible_float");
        const std::vector<float> samples = {-1.f, -0.25f, 0.f, 0.5f, 1.f};
        std::vector<uint8_t> bytes;
        bytes.reserve(samples.size() * sizeof(float));
        for (float sample : samples) {
            uint32_t bits = 0;
            std::memcpy(&bits, &sample, sizeof(bits));
            bytes.push_back(uint8_t(bits));
            bytes.push_back(uint8_t(bits >> 8));
            bytes.push_back(uint8_t(bits >> 16));
            bytes.push_back(uint8_t(bits >> 24));
        }
        writeExtensibleFixture(path, 3, bytes, 32);

        wm::WaveFileReader reader(path.string());
        CHECK(reader.getAudioFormat() == 3);
        CHECK(reader.getSampleBitDepth() == 32);
        REQUIRE(reader.getNumFrames() == samples.size());

        std::vector<float> decoded(samples.size());
        REQUIRE(reader.readFrames(0, uint32_t(samples.size()), decoded.data()) == samples.size());
        for (size_t i = 0; i < samples.size(); ++i) {
            CHECK(decoded[i] == samples[i]);
        }

        reader.close();
        std::filesystem::remove(path);
    }
}

TEST_CASE("WAVE_FORMAT_EXTENSIBLE rejects unsupported sub-formats")
{
    const auto path = tempWavPath("extensible_unsupported");
    writeExtensibleFixture(path, 2, {0, 0}, 16);

    CHECK_THROWS(wm::WaveFileReader(path.string()));

    std::filesystem::remove(path);
}

// ============================================================================
// Random-access reads: unaligned mid-range and single-frame reads
// ============================================================================

TEST_CASE("unaligned mid-range readFrames() matches the equivalent oracle slice")
{
    wm::Wave wave = makeTestWave(1000, 2, 16, 48000, true);
    const auto path = tempWavPath("midrange");
    wave.saveAs(path.string());

    wm::Wave oracle(path.string());
    wm::WaveFileReader reader(path.string());

    const uint16_t channels = oracle.getNumChannels();
    const uint32_t startFrame = 137;
    const uint32_t frameCount = 253;
    std::vector<float> buffer(size_t(frameCount) * channels);

    const uint32_t got = reader.readFrames(startFrame, frameCount, buffer.data());
    REQUIRE(got == frameCount);

    const float* expected = oracle.constAudioData() + size_t(startFrame) * channels;
    for (size_t i = 0; i < buffer.size(); ++i) {
        CHECK(buffer[i] == expected[i]);
    }

    reader.close(); // Windows can't delete a file with an open handle.
    std::filesystem::remove(path);
}

TEST_CASE("independent readers over the same file do not share a seek position")
{
    wm::Wave wave = makeTestWave(200, 1, 16, 44100, true);
    const auto path = tempWavPath("independent");
    wave.saveAs(path.string());

    wm::Wave oracle(path.string());
    wm::WaveFileReader readerA(path.string());
    wm::WaveFileReader readerB(path.string());

    std::vector<float> bufferA(10);
    std::vector<float> bufferB(10);

    // Interleave calls on the two readers to prove neither perturbs the other's position.
    REQUIRE(readerB.readFrames(150, 10, bufferB.data()) == 10);
    REQUIRE(readerA.readFrames(0, 10, bufferA.data()) == 10);
    REQUIRE(readerB.readFrames(50, 10, bufferB.data()) == 10);

    const float* expectedAtZero = oracle.constAudioData();
    const float* expectedAt50 = oracle.constAudioData() + 50;
    for (int i = 0; i < 10; ++i) {
        CHECK(bufferA[i] == expectedAtZero[i]);
        CHECK(bufferB[i] == expectedAt50[i]);
    }

    readerA.close(); // Windows can't delete a file with an open handle.
    readerB.close();
    std::filesystem::remove(path);
}

// ============================================================================
// Boundary reads
// ============================================================================

TEST_CASE("readFrames() boundary conditions")
{
    wm::Wave wave = makeTestWave(100, 2, 16, 44100, true);
    const auto path = tempWavPath("boundary");
    wave.saveAs(path.string());

    wm::WaveFileReader reader(path.string());
    const uint16_t channels = reader.getNumChannels();
    const uint32_t numFrames = reader.getNumFrames();
    REQUIRE(numFrames == 100);

    SUBCASE("frameCount of zero reads nothing")
    {
        std::vector<float> buffer(channels, -2.f);
        CHECK(reader.readFrames(0, 0, buffer.data()) == 0);
    }

    SUBCASE("startFrame exactly at frame count reads nothing")
    {
        std::vector<float> buffer(channels, -2.f);
        CHECK(reader.readFrames(numFrames, 10, buffer.data()) == 0);
    }

    SUBCASE("startFrame past frame count reads nothing")
    {
        std::vector<float> buffer(channels, -2.f);
        CHECK(reader.readFrames(numFrames + 1000, 10, buffer.data()) == 0);
    }

    SUBCASE("a read past the end is truncated, not zero-filled")
    {
        std::vector<float> buffer(size_t(channels) * 10, -2.f);
        const uint32_t got = reader.readFrames(numFrames - 3, 10, buffer.data());
        CHECK(got == 3);
    }

    SUBCASE("reading exactly the last frame succeeds")
    {
        wm::Wave oracle(path.string());
        std::vector<float> buffer(channels);
        REQUIRE(reader.readFrames(numFrames - 1, 1, buffer.data()) == 1);
        const float* expected = oracle.constAudioData() + size_t(numFrames - 1) * channels;
        for (uint16_t ch = 0; ch < channels; ++ch) {
            CHECK(buffer[ch] == expected[ch]);
        }
    }

    reader.close(); // Windows can't delete a file with an open handle.
    std::filesystem::remove(path);
}

// ============================================================================
// Odd chunk layout: extra chunks before/between fmt and data, with padding
// ============================================================================

// Note: Wave and WaveFileReader share wm::detail::parseWaveFile() (see
// WaveManipp/src/WaveFileFormat.hpp), so comparing this handcrafted fixture's decode
// against a freshly-opened wm::Wave no longer proves header-parsing correctness
// independently -- a shared parser bug would pass both. The assertions below check
// against explicit, independently-known expected values instead; the wm::Wave
// comparison is kept only as a secondary integration sanity check.
TEST_CASE("an odd chunk layout with extra chunks and padding parses correctly")
{
    // Hand-assembled RIFF/WAVE file: an odd-sized "JUNK" chunk (exercises the
    // padding-byte skip) before "fmt ", then a "LIST" chunk between "fmt " and
    // "data" -- both should be skipped by the generic chunk walk.
    const uint16_t channels = 1;
    const uint32_t sampleRate = 44100;
    const uint16_t bitsPerSample = 16;
    const std::vector<int16_t> samples = {0, 1000, -1000, 32000, -32000, 12345, -12345, 500};
    const uint32_t dataSize = uint32_t(samples.size() * sizeof(int16_t));

    const std::vector<uint8_t> junkPayload = {1, 2, 3}; // odd size -> one padding byte follows
    const std::vector<uint8_t> listPayload = {'I', 'N', 'F', 'O', 'x', 'y'}; // even size, no padding

    const uint32_t fmtChunkSize = 16;
    const uint32_t riffChunkSize = 4 /*"WAVE"*/
        + 8 + uint32_t(junkPayload.size()) + uint32_t(junkPayload.size() % 2)
        + 8 + fmtChunkSize
        + 8 + uint32_t(listPayload.size()) + uint32_t(listPayload.size() % 2)
        + 8 + dataSize;

    const auto path = tempWavPath("oddchunks");
    {
        std::ofstream out(path, std::ios::binary);
        auto writeTag = [&](const char* tag) { out.write(tag, 4); };
        auto writeU32 = [&](uint32_t v) { out.write(reinterpret_cast<const char*>(&v), 4); };
        auto writeU16 = [&](uint16_t v) { out.write(reinterpret_cast<const char*>(&v), 2); };

        writeTag("RIFF");
        writeU32(riffChunkSize);
        writeTag("WAVE");

        writeTag("JUNK");
        writeU32(uint32_t(junkPayload.size()));
        out.write(reinterpret_cast<const char*>(junkPayload.data()), junkPayload.size());
        out.put('\0'); // padding byte for the odd-sized chunk

        writeTag("fmt ");
        writeU32(fmtChunkSize);
        writeU16(1); // PCM
        writeU16(channels);
        writeU32(sampleRate);
        writeU32(sampleRate * channels * bitsPerSample / 8); // byte rate
        writeU16(channels * bitsPerSample / 8); // block align
        writeU16(bitsPerSample);

        writeTag("LIST");
        writeU32(uint32_t(listPayload.size()));
        out.write(reinterpret_cast<const char*>(listPayload.data()), listPayload.size());

        writeTag("data");
        writeU32(dataSize);
        out.write(reinterpret_cast<const char*>(samples.data()), dataSize);
    }

    wm::WaveFileReader reader(path.string());

    // Explicit expected metadata, independent of any oracle.
    CHECK(reader.getAudioFormat() == 1); // PCM
    CHECK(reader.getNumChannels() == channels);
    CHECK(reader.getSampleRate() == sampleRate);
    CHECK(reader.getSampleBitDepth() == bitsPerSample);
    CHECK(reader.isLittleEndian());
    REQUIRE(reader.getNumFrames() == samples.size());

    std::vector<float> buffer(samples.size());
    REQUIRE(reader.readFrames(0, uint32_t(samples.size()), buffer.data()) == samples.size());

    // Explicit expected sample values: PCM16 -> normalized float is `raw / 32768.f`,
    // computed here directly from the known `samples` array, not read back from any
    // other decoder.
    for (size_t i = 0; i < samples.size(); ++i) {
        const float expected = float(samples[i]) / 32768.f;
        CHECK(buffer[i] == expected);
    }

    // Secondary integration sanity check, not an independent proof (see note above).
    wm::Wave oracle(path.string());
    checkMetadataMatches(oracle, reader);
    const float* oracleData = oracle.constAudioData();
    for (size_t i = 0; i < samples.size(); ++i) {
        CHECK(buffer[i] == oracleData[i]);
    }

    reader.close(); // Windows can't delete a file with an open handle.
    std::filesystem::remove(path);
}

// ============================================================================
// Truncated / corrupt input
// ============================================================================

TEST_CASE("a file missing the fmt and data chunks throws on open, matching wm::Wave")
{
    const auto path = tempWavPath("truncated");
    {
        std::ofstream out(path, std::ios::binary);
        out.write("RIFF", 4);
        const uint32_t riffSize = 4;
        out.write(reinterpret_cast<const char*>(&riffSize), 4);
        out.write("WAVE", 4);
        // No fmt or data chunk follows.
    }

    CHECK_THROWS(wm::Wave(path.string()));
    CHECK_THROWS(wm::WaveFileReader(path.string()));

    std::filesystem::remove(path);
}

TEST_CASE("a data chunk whose declared payload exceeds the file is rejected")
{
    const auto path = tempWavPath("truncated_payload");
    {
        std::ofstream out(path, std::ios::binary);
        out.write("RIFF", 4);
        writeU32(out, 38, true); // Physical RIFF length: 46 bytes.
        out.write("WAVEfmt ", 8);
        writeU32(out, 16, true);
        writeU16(out, 1, true);
        writeU16(out, 1, true);
        writeU32(out, 44100, true);
        writeU32(out, 44100 * 2, true);
        writeU16(out, 2, true);
        writeU16(out, 16, true);
        out.write("data", 4);
        writeU32(out, 4, true); // Claims two samples...
        writeU16(out, 1234, true); // ...but only one is present.
    }

    CHECK_THROWS(wm::Wave(path.string()));
    CHECK_THROWS(wm::WaveFileReader(path.string()));

    std::filesystem::remove(path);
}

TEST_CASE("a file that is not RIFF/RIFX throws on open")
{
    const auto path = tempWavPath("notriff");
    {
        std::ofstream out(path, std::ios::binary);
        out.write("FOOB", 4);
        const uint32_t size = 4;
        out.write(reinterpret_cast<const char*>(&size), 4);
        out.write("WAVE", 4);
    }

    CHECK_THROWS(wm::WaveFileReader(path.string()));

    std::filesystem::remove(path);
}

TEST_CASE("readFrames() on an empty reader (never opened) reads nothing")
{
    wm::WaveFileReader reader;
    CHECK_FALSE(reader.isOpen());
    std::vector<float> buffer(10);
    CHECK(reader.readFrames(0, 10, buffer.data()) == 0);
}

TEST_CASE("move construction and move assignment transfer the open file")
{
    wm::Wave wave = makeTestWave(50, 1, 16, 44100, true);
    const auto path = tempWavPath("move");
    wave.saveAs(path.string());

    wm::WaveFileReader original(path.string());
    wm::WaveFileReader moved(std::move(original));
    CHECK_FALSE(original.isOpen());
    CHECK(moved.isOpen());
    CHECK(moved.getNumFrames() == 50);

    std::vector<float> buffer(1);
    CHECK(moved.readFrames(0, 1, buffer.data()) == 1);

    wm::WaveFileReader assignedTo;
    assignedTo = std::move(moved);
    CHECK_FALSE(moved.isOpen());
    CHECK(assignedTo.isOpen());
    CHECK(assignedTo.readFrames(0, 1, buffer.data()) == 1);

    assignedTo.close(); // Windows can't delete a file with an open handle.
    std::filesystem::remove(path);
}
