#pragma once

/*!
 * \file Fft.hpp
 * \brief FFT plans, transforms, convolution, and Wave spectrum helpers.
 */

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

/*!
 * \brief Returns true when \p n is a nonzero power of two.
 * \param n Value to test.
 */
bool isPowerOfTwo(std::size_t n);
/*!
 * \brief Returns true when \p n is a nonzero power of four.
 * \param n Value to test.
 */
bool isPowerOfFour(std::size_t n);

/*!
 * \brief Reusable FFT plan with precomputed twiddle tables.
 *
 * \details Plan is useful when transforming many vectors of the same size. The
 * transform operates in place on separate real and imaginary arrays with matching
 * lengths equal to size().
 */
class Plan {
public:
    /*!
     * \brief Creates a plan for a given transform size and algorithm.
     * \param size      Transform size in samples. Radix2/Radix4 require a matching power.
     * \param algorithm Which FFT algorithm to use; Auto selects the best fit for \p size.
     */
    explicit Plan(std::size_t size, Algorithm algorithm = Algorithm::Auto);

    /*! \brief Returns the planned transform size. */
    std::size_t size() const;
    /*! \brief Returns the algorithm selected by the constructor. */
    Algorithm algorithm() const;

    /*!
     * \brief Performs an in-place forward complex FFT.
     * \param real In/out array of real parts; length must equal size().
     * \param imag In/out array of imaginary parts; length must equal size().
     */
    void transform(std::vector<double>& real, std::vector<double>& imag) const;
    /*!
     * \brief Performs an in-place inverse complex FFT.
     * \param real  In/out array of real parts; length must equal size().
     * \param imag  In/out array of imaginary parts; length must equal size().
     * \param scale When true, divides each output element by size() after the transform.
     */
    void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                          bool scale = true) const;

private:
    std::size_t m_size;
    Algorithm m_algorithm;
    unsigned int m_levels;
    std::vector<double> m_cosTable;
    std::vector<double> m_sinTable;
};

/*!
 * \brief Performs an in-place forward complex FFT.
 * \param real      In/out real parts. Length determines the transform size.
 * \param imag      In/out imaginary parts. Must match the length of \p real.
 * \param algorithm Which FFT algorithm to use; Auto selects based on size.
 */
void transform(std::vector<double>& real, std::vector<double>& imag,
               Algorithm algorithm = Algorithm::Auto);
/*!
 * \brief Performs an in-place inverse complex FFT.
 * \param real      In/out real parts. Length determines the transform size.
 * \param imag      In/out imaginary parts. Must match the length of \p real.
 * \param scale     When true, divides each output element by the transform size.
 * \param algorithm Which FFT algorithm to use; Auto selects based on size.
 */
void inverseTransform(std::vector<double>& real, std::vector<double>& imag,
                      bool scale = true, Algorithm algorithm = Algorithm::Auto);

/*!
 * \brief Performs an in-place radix-2 FFT.
 * \param real In/out real parts; length must be a power of two.
 * \param imag In/out imaginary parts; must match the length of \p real.
 */
void transformRadix2(std::vector<double>& real, std::vector<double>& imag);
/*!
 * \brief Performs an in-place radix-4 FFT.
 * \param real In/out real parts; length must be a power of four.
 * \param imag In/out imaginary parts; must match the length of \p real.
 */
void transformRadix4(std::vector<double>& real, std::vector<double>& imag);
/*!
 * \brief Performs an in-place Bluestein FFT for arbitrary input sizes.
 * \param real In/out real parts.
 * \param imag In/out imaginary parts; must match the length of \p real.
 */
void transformBluestein(std::vector<double>& real, std::vector<double>& imag);

/*!
 * \brief Computes real circular convolution of two vectors.
 * \param x   First real-valued input sequence.
 * \param y   Second real-valued input sequence; must have the same length as \p x.
 * \param out Output sequence; resized to match the input length.
 */
void convolve(const std::vector<double>& x, const std::vector<double>& y,
              std::vector<double>& out);
/*!
 * \brief Computes complex circular convolution of two complex sequences.
 * \param xreal   Real part of the first input; determines the transform size.
 * \param ximag   Imaginary part of the first input; must match length of \p xreal.
 * \param yreal   Real part of the second input; must match length of \p xreal.
 * \param yimag   Imaginary part of the second input; must match length of \p xreal.
 * \param outreal Real part of the output; resized to match the input length.
 * \param outimag Imaginary part of the output; resized to match the input length.
 */
void convolve(const std::vector<double>& xreal, const std::vector<double>& ximag,
              const std::vector<double>& yreal, const std::vector<double>& yimag,
              std::vector<double>& outreal, std::vector<double>& outimag);

/*!
 * \brief Extracts a Wave channel into real FFT input with a zero-filled imaginary array.
 * \param wave            Source wave.
 * \param real            Output real array filled with channel samples.
 * \param imag            Output imaginary array filled with zeros, same length as \p real.
 * \param channel         Zero-based channel index to extract.
 * \param startFrame      First frame to include (0 = beginning of wave).
 * \param frameCount      Number of frames to extract; 0 means all frames from \p startFrame.
 * \param padToPowerOfTwo When true, zero-pads the output arrays to the next power of two.
 */
void extractChannel(const Wave& wave, std::vector<double>& real, std::vector<double>& imag,
                    int channel = 0, uint32_t startFrame = 0, uint32_t frameCount = 0,
                    bool padToPowerOfTwo = true);

/*!
 * \brief Returns a single-sided linear magnitude spectrum for one Wave channel.
 * \param wave            Source wave to analyze.
 * \param channel         Zero-based channel index.
 * \param startFrame      First frame to include (0 = beginning of wave).
 * \param frameCount      Number of frames to transform; 0 means all frames from \p startFrame.
 * \param padToPowerOfTwo When true, zero-pads input to the next power of two before transforming.
 * \param algorithm       Which FFT algorithm to use; Auto selects based on the padded size.
 * \return Single-sided magnitude spectrum with (N/2 + 1) bins, where N is the padded frame count.
 */
std::vector<double> magnitudeSpectrum(const Wave& wave, int channel = 0,
                                      uint32_t startFrame = 0, uint32_t frameCount = 0,
                                      bool padToPowerOfTwo = true,
                                      Algorithm algorithm = Algorithm::Auto);

} // namespace wm::dsp::fft

#endif // !FFT_H
