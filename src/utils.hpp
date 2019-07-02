#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

[[noreturn]] inline void throwError(const std::string& errorMsg, const std::string& funcName = "")
{
    std::string prefix = funcName != "" ? funcName + ": " : "";
    std::cerr << prefix << errorMsg << std::endl;
    throw errorMsg;
}

/**
 * @brief WaveManip
 *
 * \ingroup wm
 */
namespace wm {


}

#endif // !UTILS_H
