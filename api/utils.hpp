#pragma once

#ifndef UTILS_H
#define UTILS_H

#ifdef WAVEMANIPP_DLL
#define WAVEMANIPPAPI __declspec(dllexport)
#else
#define WAVEMANIPPAPI
#endif // WAVEMANIPP_DLL

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef __GNUC__
#include <cstring>
#include <cmath>
#endif

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
