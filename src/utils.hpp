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

    /*!
     * \brief Checks whether the current system architecture is big-endian.
     * \result <b>true</b> if the CPU is big-endian, <b>false</b> if it is little-endian
     */


    /*!
     * \brief Interpets the char array.
     */
    //template <typename T>
    //T fromLittleEndian(const char* src)
    //{
    //    if (!std::is_integral<T>::value) {
    //        throwError("Wrong Type");
    //    }
    //    size_t count = sizeof(T);
    //    T res = 0;
    //    if (isBigEndian()) {
    //        char* pRes = reinterpret_cast<char*>(&res);
    //        std::reverse_copy(src, src + count, pRes);
    //    } else {
    //        std::memcpy(&res, src, count);
    //    }
    //    return res;
    //}

    //template <typename T>
    //void fromLittleEndian(T& val)
    //{
    //    if (!std::is_integral<T>::value) {
    //        throwError("Wrong Type");
    //    }
    //    size_t count = sizeof(T);
    //    if (true) {
    //        T res = 0;
    //        char* pVal = reinterpret_cast<char*>(&val);
    //        std::reverse(pVal, pVal + count);
    //    }
    //}

}

#endif // !UTILS_H
