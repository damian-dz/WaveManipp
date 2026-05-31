#pragma once

#ifndef FFT_H
#define FFT_H

#include "Wave.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

/*!
 * \namespace wm::dsp::fft
 * \brief FFT, convolution, and Wave spectrum helpers.
 */
namespace wm::dsp::fft {

/*! \brief FFT implementation choice. */
enum class Algorithm {
    /*! \brief Choose Radix4, Radix2, or Bluestein based on transform size. */
    Auto,
    /*! \brief Cooley-Tukey radix-2 FFT for power-of-two sizes. */
    Radix2,
    /*! \brief Cooley-Tukey radix-4 FFT for power-of-four sizes. */
    Radix4,
    /*! \brief Bluestein FFT for arbitrary sizes. */
    Bluestein
};

/*! \brief Returns true when n is a nonzero power of two. */
bool isPowerOfTwo(std::size_t n);
/*! \brief Returns true when n is a nonzero power of four. */
bool isPowerOfFour(std::size_t n);

/*!
 * \brief Reusable FFT plan with precomputed tables.
 *
 * \details Plan is useful when transforming many vectors of the same size. The
 * transform operates in place on separate real and imaginary arrays with matching
 * lengths equal to size().
 */
class Plan {
public:
    /*! \brief Creates a plan for a transform size and algorithm. */
    explicit Plan(std::size_t size, Algorithm algorithm = Algorithm::Auto);

    /*! \brief Returns the planned transform size. */
    std::size_t size() const;
    /*! \brief Returns the algorithm selected by the constructor. */
    Algorithm algorithm() const;

    /*! \brief Performs an in-place forward complex FFT. */
    void transform(std::vector<double>& real, std::vector<double>& imag) const;
    /*! \brief Performs an in-place inverse complex FFT. */
    void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                          bool scale = true) const;

private:
    std::size_t m_size;
    Algorithm m_algorithm;
    unsigned int m_levels;
    std::vector<double> m_cosTable;
    std::vector<double> m_sinTable;
};

/*! \brief Performs an in-place forward complex FFT. */
void transform(std::vector<double>& real, std::vector<double>& imag,
               Algorithm algorithm = Algorithm::Auto);
/*! \brief Performs an in-place inverse complex FFT. */
void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                      bool scale = true, Algorithm algorithm = Algorithm::Auto);

/*! \brief Performs an in-place radix-2 FFT; input size must be a power of two. */
void transformRadix2(std::vector<double>& real, std::vector<double>& imag);
/*! \brief Performs an in-place radix-4 FFT; input size must be a power of four. */
void transformRadix4(std::vector<double>& real, std::vector<double>& imag);
/*! \brief Performs an in-place Bluestein FFT for arbitrary input sizes. */
void transformBluestein(std::vector<double>& real, std::vector<double>& imag);

/*! \brief Computes real circular convolution for vectors of matching length. */
void convolve(const std::vector<double>& x, const std::vector<double>& y,
              std::vector<double>& out);
/*! \brief Computes complex circular convolution for vectors of matching length. */
void convolve(const std::vector<double>& xreal, const std::vector<double>& ximag,
              const std::vector<double>& yreal, const std::vector<double>& yimag,
              std::vector<double>& outreal, std::vector<double>& outimag);

/*! \brief Extracts a Wave channel into real FFT input and zero-filled imaginary input. */
void extractChannel(const Wave& wave, std::vector<double>& real, std::vector<double>& imag,
                    int channel = 0, uint32_t startFrame = 0, uint32_t frameCount = 0,
                    bool padToPowerOfTwo = true);

/*! \brief Returns a single-sided linear magnitude spectrum for one Wave channel. */
std::vector<double> magnitudeSpectrum(const Wave& wave, int channel = 0,
                                      uint32_t startFrame = 0, uint32_t frameCount = 0,
                                      bool padToPowerOfTwo = true,
                                      Algorithm algorithm = Algorithm::Auto);

} // namespace wm::dsp::fft

#endif // !FFT_H
