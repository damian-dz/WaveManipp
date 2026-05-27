#include "Fft.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

std::size_t reverseBase4Digits(std::size_t x, unsigned int digits)
{
    std::size_t result = 0;
    for (unsigned int i = 0; i < digits; ++i) {
        result = result * 4 + x % 4;
        x /= 4;
    }
    return result;
}

std::size_t nextPowerOfTwo(std::size_t n)
{
    if (n == 0)
        return 0;
    std::size_t result = 1;
    while (result < n) {
        if (result > std::numeric_limits<std::size_t>::max() / 2)
            throwError("Vector too large.", "wm::dsp::fft::nextPowerOfTwo");
        result *= 2;
    }
    return result;
}

void checkSameLength(const std::vector<double>& real, const std::vector<double>& imag,
                     const std::string& funcName)
{
    if (real.size() != imag.size())
        throwError("Mismatched lengths.", funcName);
}

unsigned int radix2Levels(std::size_t n, const std::string& funcName)
{
    unsigned int levels = 0;
    for (std::size_t temp = n; temp > 1; temp >>= 1)
        ++levels;
    if (n != 0 && (std::size_t{1} << levels) != n)
        throwError("Length is not a power of 2.", funcName);
    return levels;
}

unsigned int radix4Levels(std::size_t n, const std::string& funcName)
{
    unsigned int levels = 0;
    for (std::size_t temp = n; temp > 1; temp /= 4) {
        if (temp % 4 != 0)
            throwError("Length is not a power of 4.", funcName);
        ++levels;
    }
    return levels;
}

void transformRadix2WithTables(std::vector<double>& real, std::vector<double>& imag,
                               unsigned int levels,
                               const std::vector<double>& cosTable,
                               const std::vector<double>& sinTable)
{
    const std::size_t n = real.size();
    for (std::size_t i = 0, j = 0; i < n; ++i) {
        if (j > i) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }

        std::size_t bit = n >> 1;
        while ((j & bit) != 0) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
    }

    for (std::size_t size = 2; size <= n; size *= 2) {
        const std::size_t halfsize = size / 2;
        const std::size_t tablestep = n / size;
        for (std::size_t i = 0; i < n; i += size) {
            for (std::size_t j = i, k = 0; j < i + halfsize; ++j, k += tablestep) {
                const double tpre =  real[j + halfsize] * cosTable[k] + imag[j + halfsize] * sinTable[k];
                const double tpim = -real[j + halfsize] * sinTable[k] + imag[j + halfsize] * cosTable[k];
                real[j + halfsize] = real[j] - tpre;
                imag[j + halfsize] = imag[j] - tpim;
                real[j] += tpre;
                imag[j] += tpim;
            }
        }
        if (size == n)
            break;
    }

    (void)levels;
}

void transformRadix4WithTables(std::vector<double>& real, std::vector<double>& imag,
                               unsigned int levels,
                               const std::vector<double>& cosTable,
                               const std::vector<double>& sinTable)
{
    const std::size_t n = real.size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = reverseBase4Digits(i, levels);
        if (j > i) {
            std::swap(real[i], real[j]);
            std::swap(imag[i], imag[j]);
        }
    }

    for (std::size_t size = 4; size <= n; size *= 4) {
        const std::size_t quarter = size / 4;
        const std::size_t tablestep = n / size;

        for (std::size_t i = 0; i < n; i += size) {
            for (std::size_t j = 0; j < quarter; ++j) {
                const std::size_t i0 = i + j;
                const std::size_t i1 = i0 + quarter;
                const std::size_t i2 = i1 + quarter;
                const std::size_t i3 = i2 + quarter;

                const std::size_t k1 = j * tablestep;
                const std::size_t k2 = 2 * k1;
                const std::size_t k3 = 3 * k1;

                const double r0 = real[i0];
                const double im0 = imag[i0];
                const double r1 =  real[i1] * cosTable[k1] + imag[i1] * sinTable[k1];
                const double im1 = -real[i1] * sinTable[k1] + imag[i1] * cosTable[k1];
                const double r2 =  real[i2] * cosTable[k2] + imag[i2] * sinTable[k2];
                const double im2 = -real[i2] * sinTable[k2] + imag[i2] * cosTable[k2];
                const double r3 =  real[i3] * cosTable[k3] + imag[i3] * sinTable[k3];
                const double im3 = -real[i3] * sinTable[k3] + imag[i3] * cosTable[k3];

                const double a0r = r0 + r2;
                const double a0i = im0 + im2;
                const double a1r = r0 - r2;
                const double a1i = im0 - im2;
                const double a2r = r1 + r3;
                const double a2i = im1 + im3;
                const double a3r = im1 - im3;
                const double a3i = r3 - r1;

                real[i0] = a0r + a2r;
                imag[i0] = a0i + a2i;
                real[i1] = a1r + a3r;
                imag[i1] = a1i + a3i;
                real[i2] = a0r - a2r;
                imag[i2] = a0i - a2i;
                real[i3] = a1r - a3r;
                imag[i3] = a1i - a3i;
            }
        }

        if (size == n)
            break;
    }
}

} // namespace

namespace wm::dsp::fft {

bool isPowerOfTwo(std::size_t n)
{
    return n != 0 && (n & (n - 1)) == 0;
}

bool isPowerOfFour(std::size_t n)
{
    if (n == 0)
        return false;
    while (n > 1) {
        if (n % 4 != 0)
            return false;
        n /= 4;
    }
    return true;
}

Plan::Plan(std::size_t size, Algorithm algorithm) :
    m_size(size),
    m_algorithm(algorithm),
    m_levels(0)
{
    if (m_size == 0)
        return;

    if (m_algorithm == Algorithm::Auto) {
        if (isPowerOfFour(m_size))
            m_algorithm = Algorithm::Radix4;
        else if (isPowerOfTwo(m_size))
            m_algorithm = Algorithm::Radix2;
        else
            m_algorithm = Algorithm::Bluestein;
    }

    if (m_algorithm == Algorithm::Radix2) {
        m_levels = radix2Levels(m_size, "wm::dsp::fft::Plan");
        m_cosTable.resize(m_size / 2);
        m_sinTable.resize(m_size / 2);
        for (std::size_t i = 0; i < m_size / 2; ++i) {
            m_cosTable[i] = std::cos(2.0 * pi * i / m_size);
            m_sinTable[i] = std::sin(2.0 * pi * i / m_size);
        }
    } else if (m_algorithm == Algorithm::Radix4) {
        m_levels = radix4Levels(m_size, "wm::dsp::fft::Plan");
        m_cosTable.resize(m_size);
        m_sinTable.resize(m_size);
        for (std::size_t i = 0; i < m_size; ++i) {
            m_cosTable[i] = std::cos(2.0 * pi * i / m_size);
            m_sinTable[i] = std::sin(2.0 * pi * i / m_size);
        }
    }
}

std::size_t Plan::size() const
{
    return m_size;
}

Algorithm Plan::algorithm() const
{
    return m_algorithm;
}

void Plan::transform(std::vector<double>& real, std::vector<double>& imag) const
{
    checkSameLength(real, imag, "wm::dsp::fft::Plan::transform");
    if (real.size() != m_size)
        throwError("Input length does not match plan size.", "wm::dsp::fft::Plan::transform");
    if (m_size == 0)
        return;

    if (m_algorithm == Algorithm::Radix2)
        transformRadix2WithTables(real, imag, m_levels, m_cosTable, m_sinTable);
    else if (m_algorithm == Algorithm::Radix4)
        transformRadix4WithTables(real, imag, m_levels, m_cosTable, m_sinTable);
    else
        transformBluestein(real, imag);
}

void Plan::inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                            bool scale) const
{
    checkSameLength(real, imag, "wm::dsp::fft::Plan::inverseTransform");
    transform(imag, real);
    if (!scale || real.empty())
        return;

    const double factor = 1.0 / static_cast<double>(real.size());
    for (std::size_t i = 0; i < real.size(); ++i) {
        real[i] *= factor;
        imag[i] *= factor;
    }
}

void transform(std::vector<double>& real, std::vector<double>& imag, Algorithm algorithm)
{
    checkSameLength(real, imag, "wm::dsp::fft::transform");
    const std::size_t n = real.size();
    if (n == 0)
        return;

    Plan(n, algorithm).transform(real, imag);
}

void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                      bool scale, Algorithm algorithm)
{
    checkSameLength(real, imag, "wm::dsp::fft::inverseTransform");
    Plan(real.size(), algorithm).inverseTransform(real, imag, scale);
}

void transformRadix2(std::vector<double>& real, std::vector<double>& imag)
{
    checkSameLength(real, imag, "wm::dsp::fft::transformRadix2");
    const std::size_t n = real.size();
    if (n == 0)
        return;

    Plan(n, Algorithm::Radix2).transform(real, imag);
}

void transformRadix4(std::vector<double>& real, std::vector<double>& imag)
{
    checkSameLength(real, imag, "wm::dsp::fft::transformRadix4");
    const std::size_t n = real.size();
    if (n == 0)
        return;

    Plan(n, Algorithm::Radix4).transform(real, imag);
}

void transformBluestein(std::vector<double>& real, std::vector<double>& imag)
{
    checkSameLength(real, imag, "wm::dsp::fft::transformBluestein");
    const std::size_t n = real.size();
    if (n == 0)
        return;
    if (n > (std::numeric_limits<std::size_t>::max() - 1) / 2)
        throwError("Vector too large.", "wm::dsp::fft::transformBluestein");

    const std::size_t m = nextPowerOfTwo(n * 2 + 1);

    std::vector<double> cosTable(n);
    std::vector<double> sinTable(n);
    const unsigned long long modulus = static_cast<unsigned long long>(n) * 2ull;
    for (std::size_t i = 0; i < n; ++i) {
        const auto ii = static_cast<unsigned long long>(i);
        const double angle = pi * static_cast<double>((ii * ii) % modulus) / static_cast<double>(n);
        cosTable[i] = std::cos(angle);
        sinTable[i] = std::sin(angle);
    }

    std::vector<double> areal(m);
    std::vector<double> aimag(m);
    for (std::size_t i = 0; i < n; ++i) {
        areal[i] =  real[i] * cosTable[i] + imag[i] * sinTable[i];
        aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
    }

    std::vector<double> breal(m);
    std::vector<double> bimag(m);
    breal[0] = cosTable[0];
    bimag[0] = sinTable[0];
    for (std::size_t i = 1; i < n; ++i) {
        breal[i] = breal[m - i] = cosTable[i];
        bimag[i] = bimag[m - i] = sinTable[i];
    }

    std::vector<double> creal(m);
    std::vector<double> cimag(m);
    convolve(areal, aimag, breal, bimag, creal, cimag);

    for (std::size_t i = 0; i < n; ++i) {
        real[i] =  creal[i] * cosTable[i] + cimag[i] * sinTable[i];
        imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
    }
}

void convolve(const std::vector<double>& x, const std::vector<double>& y,
              std::vector<double>& out)
{
    if (x.size() != y.size() || x.size() != out.size())
        throwError("Mismatched lengths.", "wm::dsp::fft::convolve");

    std::vector<double> ximag(x.size());
    std::vector<double> yimag(y.size());
    std::vector<double> outimag(out.size());
    convolve(x, ximag, y, yimag, out, outimag);
}

void convolve(const std::vector<double>& xreal, const std::vector<double>& ximag,
              const std::vector<double>& yreal, const std::vector<double>& yimag,
              std::vector<double>& outreal, std::vector<double>& outimag)
{
    if (xreal.size() != ximag.size() || xreal.size() != yreal.size() ||
        yreal.size() != yimag.size() || xreal.size() != outreal.size() ||
        outreal.size() != outimag.size()) {
        throwError("Mismatched lengths.", "wm::dsp::fft::convolve");
    }

    std::vector<double> xr(xreal);
    std::vector<double> xi(ximag);
    std::vector<double> yr(yreal);
    std::vector<double> yi(yimag);

    transform(xr, xi);
    transform(yr, yi);
    for (std::size_t i = 0; i < xr.size(); ++i) {
        const double temp = xr[i] * yr[i] - xi[i] * yi[i];
        xi[i] = xi[i] * yr[i] + xr[i] * yi[i];
        xr[i] = temp;
    }
    inverseTransform(xr, xi);

    for (std::size_t i = 0; i < xr.size(); ++i) {
        outreal[i] = xr[i];
        outimag[i] = xi[i];
    }
}

void extractChannel(const Wave& wave, std::vector<double>& real, std::vector<double>& imag,
                    int channel, uint32_t startFrame, uint32_t frameCount,
                    bool padToPowerOfTwo)
{
    if (wave.isEmpty()) {
        real.clear();
        imag.clear();
        return;
    }
    if (channel < 0 || channel >= static_cast<int>(wave.getNumChannels()))
        throwError("Channel index out of range.", "wm::dsp::fft::extractChannel");
    if (startFrame > wave.getNumFrames())
        throwError("Start frame is out of range.", "wm::dsp::fft::extractChannel");

    const uint32_t available = wave.getNumFrames() - startFrame;
    const uint32_t count = frameCount == 0 ? available : std::min(frameCount, available);
    const std::size_t outputSize = padToPowerOfTwo ? nextPowerOfTwo(count) : count;

    real.assign(outputSize, 0.0);
    imag.assign(outputSize, 0.0);
    for (uint32_t i = 0; i < count; ++i)
        real[i] = wave(startFrame + i, channel);
}

std::vector<double> magnitudeSpectrum(const Wave& wave, int channel,
                                      uint32_t startFrame, uint32_t frameCount,
                                      bool padToPowerOfTwo, Algorithm algorithm)
{
    std::vector<double> real;
    std::vector<double> imag;
    extractChannel(wave, real, imag, channel, startFrame, frameCount, padToPowerOfTwo);
    transform(real, imag, algorithm);

    const std::size_t bins = real.empty() ? 0 : real.size() / 2 + 1;
    std::vector<double> magnitude(bins);
    const double scale = real.empty() ? 0.0 : 1.0 / static_cast<double>(real.size());
    for (std::size_t i = 0; i < bins; ++i) {
        const double value = std::sqrt(real[i] * real[i] + imag[i] * imag[i]) * scale;
        magnitude[i] = (i == 0 || (real.size() % 2 == 0 && i == bins - 1)) ? value : 2.0 * value;
    }
    return magnitude;
}

} // namespace wm::dsp::fft

/*
 * FFT and convolution routines adapted from:
 * Free FFT and convolution (C++) by Project Nayuki
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 *
 * MIT License
 */
