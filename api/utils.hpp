#pragma once

/*!
 * \file utils.hpp
 * \brief Shared constants, the SampleRate enum, error helper, and debug macros.
 */

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

/*! \def PRINT(var)  Prints \p var followed by a newline to stdout. */
#define PRINT(var) (std::cout << var << std::endl);
/*! \def PRINTN(var) Prints "var: value" followed by a newline to stdout. */
#define PRINTN(var) (std::cout << #var << ": " << var << std::endl);

/*! \brief Common audio sample rates in Hz. */
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

/*! \brief Default intermediate I/O buffer size in bytes. */
constexpr uint32_t default_buffer_size = 24576;
/*! \brief Default sample rate in Hz. */
constexpr uint32_t default_sample_rate = static_cast<uint32_t>(SampleRate::SR_44100);

/*!
 * \brief Prints an error message and throws it as a std::string exception payload.
 * \param errorMsg Human-readable description of the error.
 * \param funcName Optional caller name prepended to the printed message.
 */
[[noreturn]] inline void throwError(const std::string& errorMsg, const std::string& funcName = "")
{
    std::string prefix = funcName != "" ? funcName + ": " : "";
    std::cerr << prefix << errorMsg << std::endl;
    throw errorMsg;
}

/*! \namespace wm \brief Root namespace for the WaveManipp library. */
namespace wm {
}

#endif // !UTILS_H
