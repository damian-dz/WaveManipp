#pragma once

/*!
 * \file PcmCodec.hpp
 * \brief Private, sample-type-independent integer PCM conversion helpers.
 */

#ifndef PCM_CODEC_H
#define PCM_CODEC_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace wm::detail {

/*!
 * \brief Quantizes a normalized floating-point sample to signed integer PCM.
 * \tparam Sample Floating-point source type (`float` or `double`).
 * \tparam Bits   Signed PCM bit depth, from 2 through 32.
 * \details Uses saturating round-to-nearest with halfway values rounded away
 * from zero. NaN maps to silence; infinities and out-of-range finite values
 * saturate to the corresponding endpoint. Dither is intentionally outside
 * this low-level deterministic conversion.
 */
template <typename Sample, unsigned Bits>
int32_t quantizeSignedPcm(Sample value) noexcept
{
    static_assert(std::is_floating_point_v<Sample>, "PCM source samples must be floating point");
    static_assert(Bits >= 2 && Bits <= 32, "signed PCM bit depth must be between 2 and 32");

    constexpr int64_t scale = int64_t{1} << (Bits - 1);
    constexpr int64_t minimum = -scale;
    constexpr int64_t maximum = scale - 1;

    if (std::isnan(value)) {
        return 0;
    }
    if (value <= Sample{-1}) {
        return static_cast<int32_t>(minimum);
    }
    if (value >= Sample{1}) {
        return static_cast<int32_t>(maximum);
    }

    const Sample scaled = value * static_cast<Sample>(scale);
    const int64_t rounded = static_cast<int64_t>(std::round(scaled));
    return static_cast<int32_t>(std::clamp(rounded, minimum, maximum));
}

/*!
 * \brief Converts a signed integer PCM code to a normalized floating-point sample.
 * \tparam Sample Floating-point destination type (`float` or `double`).
 * \tparam Bits   Signed PCM bit depth, from 2 through 32.
 */
template <typename Sample, unsigned Bits>
Sample dequantizeSignedPcm(int32_t value) noexcept
{
    static_assert(std::is_floating_point_v<Sample>, "PCM destination samples must be floating point");
    static_assert(Bits >= 2 && Bits <= 32, "signed PCM bit depth must be between 2 and 32");

    constexpr int64_t scale = int64_t{1} << (Bits - 1);
    constexpr int64_t minimum = -scale;
    constexpr int64_t maximum = scale - 1;
    const int64_t clamped = std::clamp<int64_t>(value, minimum, maximum);
    return static_cast<Sample>(clamped) / static_cast<Sample>(scale);
}

//! Encodes one normalized sample as unsigned 8-bit WAVE PCM.
template <typename Sample>
uint8_t encodePcm8Sample(Sample value) noexcept
{
    return static_cast<uint8_t>(quantizeSignedPcm<Sample, 8>(value) + 128);
}

//! Decodes one unsigned 8-bit WAVE PCM code to a normalized sample.
template <typename Sample>
Sample decodePcm8Sample(uint8_t value) noexcept
{
    return dequantizeSignedPcm<Sample, 8>(int32_t(value) - 128);
}

//! Decodes one signed 24-bit PCM sample stored in file byte order.
template <typename Sample>
Sample decodePcm24Sample(const uint8_t* bytes, bool isLittleEndian) noexcept
{
    const uint32_t value = isLittleEndian ?
        uint32_t(bytes[0]) | (uint32_t(bytes[1]) << 8) | (uint32_t(bytes[2]) << 16) :
        uint32_t(bytes[2]) | (uint32_t(bytes[1]) << 8) | (uint32_t(bytes[0]) << 16);
    const int32_t signedValue = (value & 0x00800000u) != 0 ?
        static_cast<int32_t>(value) - 0x01000000 : static_cast<int32_t>(value);
    return dequantizeSignedPcm<Sample, 24>(signedValue);
}

//! Encodes one normalized sample as signed 24-bit PCM in file byte order.
template <typename Sample>
void encodePcm24Sample(Sample value, bool isLittleEndian, uint8_t* bytes) noexcept
{
    const uint32_t bits = static_cast<uint32_t>(quantizeSignedPcm<Sample, 24>(value)) & 0x00FFFFFFu;
    if (isLittleEndian) {
        bytes[0] = uint8_t(bits);
        bytes[1] = uint8_t(bits >> 8);
        bytes[2] = uint8_t(bits >> 16);
    } else {
        bytes[0] = uint8_t(bits >> 16);
        bytes[1] = uint8_t(bits >> 8);
        bytes[2] = uint8_t(bits);
    }
}

} // namespace wm::detail

#endif // !PCM_CODEC_H
