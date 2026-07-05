#include <doctest/doctest.h>

#include "PcmCodec.hpp"
#include "Wave.hpp"
#include "WaveFileReader.hpp"
#include "WaveFileWriter.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

std::filesystem::path writerTempPath(const char* suffix)
{
    return std::filesystem::temp_directory_path() / (std::string("wavemanipp_writer_test_") + suffix + ".wav");
}

uint32_t readU32(const std::filesystem::path& path, uint64_t offset, bool littleEndian)
{
    std::ifstream input(path, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(offset));
    uint8_t bytes[4]{};
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    if (input.gcount() != sizeof(bytes)) {
        throw std::runtime_error("failed to read checkpointed WAVE header");
    }
    return littleEndian ?
        uint32_t(bytes[0]) | (uint32_t(bytes[1]) << 8) | (uint32_t(bytes[2]) << 16) | (uint32_t(bytes[3]) << 24) :
        uint32_t(bytes[3]) | (uint32_t(bytes[2]) << 8) | (uint32_t(bytes[1]) << 16) | (uint32_t(bytes[0]) << 24);
}

float expectedDecoded(float input, const wm::WaveFileWriterFormat& format)
{
    if (format.encoding == wm::WaveFileEncoding::IeeeFloat) {
        return input;
    }
    if (format.bitsPerSample == 8) {
        const int32_t code = wm::detail::quantizeSignedPcm<float, 8>(input);
        return wm::detail::dequantizeSignedPcm<float, 8>(code);
    }
    if (format.bitsPerSample == 16) {
        const int32_t code = wm::detail::quantizeSignedPcm<float, 16>(input);
        return wm::detail::dequantizeSignedPcm<float, 16>(code);
    }
    const int32_t code = wm::detail::quantizeSignedPcm<float, 24>(input);
    return wm::detail::dequantizeSignedPcm<float, 24>(code);
}

} // namespace

TEST_CASE("WaveFileWriter incrementally writes every supported encoding and endianness")
{
    struct FormatCase {
        wm::WaveFileEncoding encoding;
        uint16_t bitsPerSample;
    };
    const std::vector<FormatCase> formats = {
        {wm::WaveFileEncoding::PcmInteger, 8},
        {wm::WaveFileEncoding::PcmInteger, 16},
        {wm::WaveFileEncoding::PcmInteger, 24},
        {wm::WaveFileEncoding::IeeeFloat, 32}
    };
    const std::vector<float> samples = {
        -1.f, -0.75f, -0.1f, 0.f, 0.1f, 0.5f, 0.999f,
        -0.333f, 0.333f, -0.0001f, 0.0001f, 0.25f, -0.25f, 0.875f,
        -0.875f, 0.625f, -0.625f, 0.05f, -0.05f, 0.7f, -0.7f, 0.2f
    };
    constexpr uint16_t channels = 2;
    constexpr uint32_t frames = 11;

    for (const FormatCase& formatCase : formats) {
        for (bool littleEndian : {true, false}) {
            CAPTURE(formatCase.bitsPerSample);
            CAPTURE(littleEndian);
            const auto path = writerTempPath(littleEndian ? "sweep_le" : "sweep_be");
            wm::WaveFileWriterFormat format;
            format.encoding = formatCase.encoding;
            format.numChannels = channels;
            format.sampleRate = 48000;
            format.bitsPerSample = formatCase.bitsPerSample;
            format.isLittleEndian = littleEndian;

            wm::WaveFileWriter writer(path.string(), format);
            REQUIRE(writer.writeFrames(samples.data(), 3) == 3);
            REQUIRE(writer.writeFrames(samples.data() + 3 * channels, 4) == 4);
            REQUIRE(writer.writeFrames(samples.data() + 7 * channels, 4) == 4);
            CHECK(writer.getNumFramesWritten() == frames);
            writer.finalize();
            CHECK_FALSE(writer.isOpen());

            wm::WaveFileReader reader(path.string());
            CHECK(reader.getAudioFormat() == static_cast<uint16_t>(format.encoding));
            CHECK(reader.getNumChannels() == channels);
            CHECK(reader.getSampleRate() == format.sampleRate);
            CHECK(reader.getSampleBitDepth() == format.bitsPerSample);
            CHECK(reader.isLittleEndian() == littleEndian);
            REQUIRE(reader.getNumFrames() == frames);

            std::vector<float> decoded(samples.size());
            REQUIRE(reader.readFrames(0, frames, decoded.data()) == frames);
            for (size_t i = 0; i < samples.size(); ++i) {
                CHECK(decoded[i] == expectedDecoded(samples[i], format));
            }

            wm::Wave loaded(path.string());
            REQUIRE(loaded.getNumFrames() == frames);
            for (size_t i = 0; i < samples.size(); ++i) {
                CHECK(loaded.constAudioData()[i] == decoded[i]);
            }

            const uint64_t dataBytes = uint64_t(frames) * channels * format.bitsPerSample / 8;
            CHECK(std::filesystem::file_size(path) == 44 + dataBytes + (dataBytes & 1u));
            reader.close();
            std::filesystem::remove(path);
        }
    }
}

TEST_CASE("WaveFileWriter checkpoints a recoverable prefix after every block")
{
    const auto path = writerTempPath("checkpoint");
    wm::WaveFileWriterFormat format;
    format.numChannels = 1;
    format.sampleRate = 44100;
    format.bitsPerSample = 8;

    wm::WaveFileWriter writer(path.string(), format);
    const float first[] = {-1.f, 0.f, 0.5f};
    REQUIRE(writer.writeFrames(first, 3) == 3);
    CHECK(readU32(path, 4, true) == 40); // 36 + 3 data bytes + one pad byte.
    CHECK(readU32(path, 40, true) == 3);
    CHECK(std::filesystem::file_size(path) == 48);

    const float second[] = {0.75f, 0.25f};
    REQUIRE(writer.writeFrames(second, 2) == 2);
    CHECK(readU32(path, 4, true) == 42); // 36 + 5 data bytes + one pad byte.
    CHECK(readU32(path, 40, true) == 5);
    CHECK(std::filesystem::file_size(path) == 50);

    writer.close(); // Simulates stopping without an explicit finalization call.
    wm::WaveFileReader reader(path.string());
    REQUIRE(reader.getNumFrames() == 5);
    std::vector<float> decoded(5);
    REQUIRE(reader.readFrames(0, 5, decoded.data()) == 5);
    CHECK(decoded[0] == -1.f);
    CHECK(decoded[1] == 0.f);
    CHECK(decoded[2] == 0.5f);
    CHECK(decoded[3] == 0.75f);
    CHECK(decoded[4] == 0.25f);

    reader.close();
    std::filesystem::remove(path);
}

TEST_CASE("WaveFileWriter accepts double input and applies the shared PCM policy")
{
    const auto path = writerTempPath("double");
    wm::WaveFileWriterFormat format;
    format.numChannels = 1;
    format.sampleRate = 96000;
    format.bitsPerSample = 24;
    const std::vector<double> samples = {
        -2.0, -1.0, -0.50000001, 0.0, 0.50000001, 1.0, 2.0,
        std::numeric_limits<double>::quiet_NaN(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()
    };

    wm::WaveFileWriter writer(path.string(), format);
    REQUIRE(writer.writeFrames(samples.data(), uint32_t(samples.size())) == samples.size());
    writer.finalize();

    wm::WaveFileReader reader(path.string());
    std::vector<float> decoded(samples.size());
    REQUIRE(reader.readFrames(0, uint32_t(samples.size()), decoded.data()) == samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        const int32_t code = wm::detail::quantizeSignedPcm<double, 24>(samples[i]);
        const float expected = wm::detail::dequantizeSignedPcm<float, 24>(code);
        CHECK(decoded[i] == expected);
    }

    reader.close();
    std::filesystem::remove(path);
}

TEST_CASE("WaveFileWriter converts double input to 32-bit IEEE float output")
{
    const auto path = writerTempPath("double_float");
    wm::WaveFileWriterFormat format;
    format.encoding = wm::WaveFileEncoding::IeeeFloat;
    format.numChannels = 1;
    format.sampleRate = 48000;
    format.bitsPerSample = 32;
    format.isLittleEndian = false;
    const std::vector<double> samples = {-1.0, -0.1, 0.0, 1.0 / 3.0, 1.0};

    wm::WaveFileWriter writer(path.string(), format);
    REQUIRE(writer.writeFrames(samples.data(), uint32_t(samples.size())) == samples.size());
    writer.finalize();

    wm::WaveFileReader reader(path.string());
    std::vector<float> decoded(samples.size());
    REQUIRE(reader.readFrames(0, uint32_t(samples.size()), decoded.data()) == samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        CHECK(decoded[i] == static_cast<float>(samples[i]));
    }

    reader.close();
    std::filesystem::remove(path);
}

TEST_CASE("WaveFileWriter finalizes empty output and finalize is idempotent")
{
    const auto path = writerTempPath("empty");
    wm::WaveFileWriter writer(path.string(), {});
    writer.checkpoint();
    writer.finalize();
    writer.finalize();
    writer.close();

    CHECK(std::filesystem::file_size(path) == 44);
    wm::WaveFileReader reader(path.string());
    CHECK(reader.getNumFrames() == 0);
    reader.close();
    std::filesystem::remove(path);
}

TEST_CASE("WaveFileWriter move operations transfer the output stream")
{
    const auto path = writerTempPath("move");
    wm::WaveFileWriterFormat format;
    format.numChannels = 1;
    format.bitsPerSample = 16;
    const float samples[] = {-0.5f, 0.f, 0.5f};

    wm::WaveFileWriter original(path.string(), format);
    REQUIRE(original.writeFrames(samples, 1) == 1);
    wm::WaveFileWriter moved(std::move(original));
    CHECK_FALSE(original.isOpen());
    CHECK(moved.isOpen());
    REQUIRE(moved.writeFrames(samples + 1, 1) == 1);

    wm::WaveFileWriter assigned;
    assigned = std::move(moved);
    CHECK_FALSE(moved.isOpen());
    REQUIRE(assigned.writeFrames(samples + 2, 1) == 1);
    assigned.finalize();

    wm::WaveFileReader reader(path.string());
    REQUIRE(reader.getNumFrames() == 3);
    reader.close();
    std::filesystem::remove(path);
}

TEST_CASE("WaveFileWriter rejects invalid formats and invalid write state")
{
    const auto path = writerTempPath("invalid");
    wm::WaveFileWriterFormat invalid;
    invalid.numChannels = 0;
    CHECK_THROWS(wm::WaveFileWriter(path.string(), invalid));

    invalid = {};
    invalid.encoding = wm::WaveFileEncoding::PcmInteger;
    invalid.bitsPerSample = 32;
    CHECK_THROWS(wm::WaveFileWriter(path.string(), invalid));

    invalid.encoding = wm::WaveFileEncoding::IeeeFloat;
    invalid.bitsPerSample = 24;
    CHECK_THROWS(wm::WaveFileWriter(path.string(), invalid));

    wm::WaveFileWriter writer;
    const float sample = 0.f;
    CHECK_THROWS(writer.writeFrames(&sample, 1));

    wm::WaveFileWriterFormat valid;
    valid.numChannels = 1;
    writer.open(path.string(), valid);
    CHECK(writer.writeFrames(static_cast<const float*>(nullptr), 0) == 0);
    CHECK_THROWS(writer.writeFrames(static_cast<const float*>(nullptr), 1));
    writer.close();

    valid.numChannels = std::numeric_limits<uint16_t>::max();
    valid.bitsPerSample = 8;
    writer.open(path.string(), valid);
    CHECK_THROWS(writer.writeFrames(&sample, std::numeric_limits<uint32_t>::max()));
    writer.close();
    std::filesystem::remove(path);
}
