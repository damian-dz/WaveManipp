#pragma once

#ifndef FFT_H
#define FFT_H

#include "Wave.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace wm::dsp::fft {

enum class Algorithm {
    Auto,
    Radix2,
    Radix4,
    Bluestein
};

bool isPowerOfTwo(std::size_t n);
bool isPowerOfFour(std::size_t n);

class Plan {
public:
    explicit Plan(std::size_t size, Algorithm algorithm = Algorithm::Auto);

    std::size_t size() const;
    Algorithm algorithm() const;

    void transform(std::vector<double>& real, std::vector<double>& imag) const;
    void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                          bool scale = true) const;

private:
    std::size_t m_size;
    Algorithm m_algorithm;
    unsigned int m_levels;
    std::vector<double> m_cosTable;
    std::vector<double> m_sinTable;
};

void transform(std::vector<double>& real, std::vector<double>& imag,
               Algorithm algorithm = Algorithm::Auto);
void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                      bool scale = true, Algorithm algorithm = Algorithm::Auto);

void transformRadix2(std::vector<double>& real, std::vector<double>& imag);
void transformRadix4(std::vector<double>& real, std::vector<double>& imag);
void transformBluestein(std::vector<double>& real, std::vector<double>& imag);

void convolve(const std::vector<double>& x, const std::vector<double>& y,
              std::vector<double>& out);
void convolve(const std::vector<double>& xreal, const std::vector<double>& ximag,
              const std::vector<double>& yreal, const std::vector<double>& yimag,
              std::vector<double>& outreal, std::vector<double>& outimag);

void extractChannel(const Wave& wave, std::vector<double>& real, std::vector<double>& imag,
                    int channel = 0, uint32_t startFrame = 0, uint32_t frameCount = 0,
                    bool padToPowerOfTwo = true);

std::vector<double> magnitudeSpectrum(const Wave& wave, int channel = 0,
                                      uint32_t startFrame = 0, uint32_t frameCount = 0,
                                      bool padToPowerOfTwo = true,
                                      Algorithm algorithm = Algorithm::Auto);

} // namespace wm::dsp::fft

#endif // !FFT_H
