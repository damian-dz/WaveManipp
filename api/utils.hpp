#pragma once

#ifndef UTILS_H
#define UTILS_H

#ifdef WAVEMANIPP_DLL
#define WAVEMANIPPAPI __declspec(dllexport)
#else
#define WAVEMANIPPAPI
#endif // WAVEMANIPP_DLL

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef __GNUC__
#include <cstring>
#include <cmath>
#endif

#define PRINT(var) (std::cout << var << std::endl);
#define PRINTN(var) (std::cout << #var << ": " << var << std::endl);

enum class SampleRate : uint32_t {
    SR_8000 = 8000,
    SR_11025 = 11025,
    SR_16000 = 16000,
    SR_22050 = 22050,
    SR_32000 = 32000,
    SR_44100 = 44100,
    SR_48000 = 48000,
    SR_88200 = 88200,
    SR_96000 = 96000,
    SR_176400 = 176400,
    SR_192000 = 192000
};

constexpr uint32_t default_buffer_size = 24576;
constexpr uint32_t default_sample_rate = static_cast<uint32_t>(SampleRate::SR_44100);

[[noreturn]] inline void throwError(const std::string& errorMsg, const std::string& funcName = "")
{
    std::string prefix = funcName != "" ? funcName + ": " : "";
    std::cerr << prefix << errorMsg << std::endl;
    throw errorMsg;
}

/**
 * @brief WaveManipp
 *
 * \ingroup wm
 */
namespace wm {


}

#endif // !UTILS_H
