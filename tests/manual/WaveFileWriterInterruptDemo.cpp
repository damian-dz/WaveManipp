// Manual diagnostic, not built by default (see EXCLUDE_FROM_ALL in CMakeLists.txt) and
// not run by ctest: an external process kills this one's "write" mode mid-stream, then
// its "verify" mode confirms the last checkpoint left a valid, playable prefix
// (docs/AUDIO_STREAMING_PLAN.md's recoverable-prefix contract for WaveFileWriter).
//
// Qualification: this tool sleeps only between completed blocks, after each block's
// writeFrames() call (and its internal checkpoint) has already returned. Killing it
// during that sleep proves recovery from abrupt termination *between* checkpoints --
// it does not exercise termination *during* fwrite() itself (a torn write to the OS).
// That narrower case is covered separately by WaveFileWriterTests.cpp's short-write
// unit test, which simulates a failing fwrite() directly. Both matter: automated tests
// can't reliably reproduce real OS-level process termination and stdio buffering
// behavior, so this manual tool remains a useful, distinct check.
//
// Usage:
//   WaveFileWriterInterruptDemo write <path> <blockCount> <sleepMs>
//   WaveFileWriterInterruptDemo verify <path>

#include "WaveFileWriter.hpp"
#include "WaveFileReader.hpp"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

namespace {

//! Parses a positive (for blockCount) or non-negative (for sleepMs) integer argument.
//! Returns false (leaving `out` unset) on anything that isn't a clean, in-range integer.
bool parsePositiveArg(const char* text, long minValue, long& out)
{
    if (!text || text[0] == '\0') {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || errno == ERANGE || value < minValue) {
        return false;
    }
    out = value;
    return true;
}

int runWrite(const std::string& path, long blockCount, long sleepMs)
{
    wm::WaveFileWriterFormat format;
    format.encoding = wm::WaveFileEncoding::PcmInteger;
    format.numChannels = 2;
    format.sampleRate = 44100;
    format.bitsPerSample = 16;
    format.isLittleEndian = true;

    wm::WaveFileWriter writer(path, format);

    constexpr uint32_t framesPerBlock = 512;
    std::vector<float> block(size_t(framesPerBlock) * format.numChannels);
    for (long b = 0; b < blockCount; ++b) {
        for (size_t i = 0; i < block.size(); ++i) {
            const float t = float(b * framesPerBlock + i) * 0.01f;
            block[i] = 0.5f * std::sin(t);
        }
        writer.writeFrames(block.data(), framesPerBlock);
        std::fprintf(stdout, "wrote block %ld/%ld (%llu frames total)\n", b + 1, blockCount,
                     (unsigned long long)writer.getNumFramesWritten());
        std::fflush(stdout);
        std::this_thread::sleep_for(std::chrono::milliseconds(sleepMs));
    }
    writer.finalize(); // Only reached if never interrupted.
    std::fprintf(stdout, "finalized normally\n");
    return 0;
}

int runVerify(const std::string& path)
{
    try {
        wm::WaveFileReader reader(path);
        std::fprintf(stdout, "opened ok: %u frames, %u channels, %u-bit, littleEndian=%d\n",
                     reader.getNumFrames(), reader.getNumChannels(), reader.getSampleBitDepth(),
                     int(reader.isLittleEndian()));

        // Bounded-memory verification: read back in fixed-size blocks rather than
        // allocating a buffer sized to the whole file, matching WaveFileReader's own
        // bounded-read contract regardless of how large the interrupted file turns out.
        constexpr uint32_t blockFrames = 4096;
        std::vector<float> block(size_t(blockFrames) * reader.getNumChannels());
        uint64_t totalRead = 0;
        while (totalRead < reader.getNumFrames()) {
            const uint32_t wanted =
                uint32_t(std::min<uint64_t>(blockFrames, reader.getNumFrames() - totalRead));
            const uint32_t got = reader.readFrames(totalRead, wanted, block.data());
            if (got == 0) {
                break;
            }
            totalRead += got;
        }
        std::fprintf(stdout, "read back %llu/%u frames successfully\n", (unsigned long long)totalRead,
                     reader.getNumFrames());
        return (totalRead == reader.getNumFrames()) ? 0 : 2;
    } catch (const std::string& err) {
        std::fprintf(stderr, "FAILED TO OPEN: %s\n", err.c_str());
        return 1;
    } catch (...) {
        std::fprintf(stderr, "FAILED TO OPEN: unknown exception\n");
        return 1;
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s write <path> <blockCount> <sleepMs> | verify <path>\n", argv[0]);
        return 64;
    }
    const std::string mode = argv[1];
    const std::string path = argv[2];

    if (mode == "write" && argc == 5) {
        long blockCount = 0;
        long sleepMs = 0;
        if (!parsePositiveArg(argv[3], 1, blockCount)) {
            std::fprintf(stderr, "blockCount must be a positive integer, got \"%s\"\n", argv[3]);
            return 64;
        }
        if (!parsePositiveArg(argv[4], 0, sleepMs)) {
            std::fprintf(stderr, "sleepMs must be a non-negative integer, got \"%s\"\n", argv[4]);
            return 64;
        }
        return runWrite(path, blockCount, sleepMs);
    }
    if (mode == "verify" && argc == 3) {
        return runVerify(path);
    }
    std::fprintf(stderr, "invalid arguments\n");
    return 64;
}
