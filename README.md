# WaveManipp

**WaveManipp** is a C++ library for uncompressed audio data manipulation. It supports WAV files up to two channels. It does not handle audio recording or playback on its own, but is designed to cooperate with third-party playback and recording libraries.

## Table of Contents:

* [Compiling from Source](#compiling-from-source)
* [Introduction](#introduction)
* [The `Wave` Class](#the-wave-class)
  * [Supported Formats](#supported-formats)
* [The `WaveBuilder` Class](#the-wavebuilder-class)
* [The `WaveMixer` Class](#the-wavemixer-class)
* [The `wm::dsp` Namespace](#the-wmdsp-namespace)
  * [FFT Helpers](#fft-helpers)

## Compiling from Source

WaveManipp uses CMake (3.20+). To build as a static library:

```
cmake -S . -B build
cmake --build build
```

To use it as a subdirectory in another CMake project:

```cmake
add_subdirectory(WaveManipp)
target_link_libraries(your_target PRIVATE WaveManipp)
```

The public headers are in `api/`. Include the umbrella header to get everything:

```cpp
#include <WaveManipp.h>
```

Building has been tested with MSVC (Visual Studio 2022) and GCC g++ using `-std=c++17`.

#### [Back to Table of Contents](#table-of-contents)

* * *
## Introduction

The easiest way to load an uncompressed audio file into memory is to pass the path to the constructor:

```c++
wm::Wave wave("D:/AudioSamples/test.wav");
```

Alternatively, construct empty and call `open()`:

```c++
wm::Wave wave;
wave.open("D:/AudioSamples/test.wav");
```

The overloaded `std::ostream` operator prints the basic header information:

```c++
std::cout << wave << std::endl;
```

```
ChunkID: RIFF
ChunkSize: 351966
Format: WAVE
----------
Subchunk1ID: fmt
Subchunk1Size: 16
AudioFormat: 1
NumChannels: 2
SampleRate: 44100
ByteRate: 176400
BlockAlign: 4
BitsPerSample: 16
----------
Subchunk2ID: data
Subchunk2Size: 351812
```

For example, generate a sine wave and append it to the existing object:

```c++
auto sineWave = wm::Wave::generateSine(220.f, 0.f, 44100, 44100);
wave.append(sineWave);
```

Save the result with `saveAs()`:

```c++
wave.saveAs("result.wav");
```

#### [Back to Table of Contents](#table-of-contents)

* * *
## The `Wave` Class

The `Wave` class represents an audio file stored in memory. All audio data is stored as single-precision floating-point values. This uses more memory than the original bit depth in most cases, but avoids per-sample conversion on every operation.

Internally, the audio buffer is reference-counted (`std::shared_ptr<float[]>`). Copying a `Wave` is O(1) — both objects share the same buffer. A private copy is made automatically (copy-on-write) the first time a mutating method is called on a shared object, so callers do not need to manage this explicitly.

```cpp
wm::Wave a("audio.wav");
wm::Wave b = a;           // O(1), shared buffer
b.changeVolume(0.5f);     // detaches; a is unaffected
```

Raw read access is available via `constAudioData()`. Mutable access via `audioData()` triggers a detach first. Note that holding a raw pointer returned by `audioData()` across a subsequent copy of the `Wave` is undefined behaviour — the copy will share the same buffer that the pointer points into.

### Supported Formats

WaveManipp supports uncompressed RIFF/RIFX mono and stereo WAV files. Supported sample bit depths are 8, 16, 24, and 32.

#### [Back to Table of Contents](#table-of-contents)

* * *
## The `WaveBuilder` Class

The `WaveBuilder` class is useful when merging a large number of `Wave` objects sequentially. Rather than repeatedly reallocating a growing buffer, it pre-calculates the total size and allocates once when `toWave()` is called.

`WaveBuilder` holds pointers to the constituent `Wave` objects, so those objects must remain alive until `toWave()` returns.

#### [Back to Table of Contents](#table-of-contents)

* * *
## The `WaveMixer` Class

This class emulates a mixing board: multiple `Wave` objects can be combined non-sequentially across named tracks. Like `WaveBuilder`, it holds pointers to existing `Wave` objects, so none of them should be destroyed before `toWave()` is called.

A `WaveMixer` starts with no tracks. Add a track, then insert audio chunks at a given frame offset:

```c++
wm::WaveMixer waveMixer;
waveMixer.addTrack();
waveMixer.insertChunk(0, 10000, sineWave);
```

#### [Back to Table of Contents](#table-of-contents)

* * *
## The `wm::dsp` Namespace

Audio processing algorithms that operate on `Wave` objects are collected in the `wm::dsp` namespace (`Dsp.hpp` / `Dsp.cpp`). These are free functions rather than member functions of `Wave`, which keeps the `Wave` class as a data container and makes the DSP operations composable independently of the file format.

Include `<WaveManipp.h>` to get the full namespace, or `<Dsp.hpp>` directly.

```cpp
#include <WaveManipp.h>

wm::Wave wave("audio.wav");

// Queries
float p = wm::dsp::peak(wave);
float pr = wm::dsp::peakRange(wave, 0, 44100);

// Gain
wm::dsp::applyGain(wave, 0.5f);                         // halve the level
wm::dsp::applyGainRange(wave, 0, 44100, 2.0f);          // double first second only

// Normalize
wm::dsp::normalize(wave);                               // normalize to 0 dBFS
wm::dsp::normalize(wave, -3.f);                         // normalize to -3 dBFS
wm::dsp::normalizeRange(wave, 0, 44100, -6.f);          // normalize first second to -6 dBFS

// Silence, invert
wm::dsp::silence(wave, 22050, 44100);                   // silence second half of first second
wm::dsp::invertRange(wave, 0, 44100);                   // phase flip first second

// Fades (logarithmic by default; pass FadeCurve::Linear for a straight ramp)
wm::dsp::fadeIn(wave, 0, 44100);
wm::dsp::fadeOut(wave, wave.getNumFrames() - 44100, wave.getNumFrames());

// Reverse
wm::dsp::reverse(wave);                                 // reverse all channels
wm::dsp::reverseChannel(wave, 0);                       // reverse left channel only

// STFT data for spectrogram drawing
wm::dsp::StftConfig stft;
stft.windowSize = 2048;
stft.hopSize = 512;
stft.channel = 0;
auto magnitudeDb = wm::dsp::stftMagnitudeDb(wave, stft); // [time frame][frequency bin]
```

All range functions clamp `endFrame` to `wave.getNumFrames()` and are no-ops if the wave is empty or `startFrame >= endFrame`.

### FFT Helpers

WaveManipp also includes vector-based FFT and convolution utilities in `wm::dsp::fft` (`Fft.hpp` / `Fft.cpp`). The automatic transform path chooses Radix-4 for power-of-4 lengths, Radix-2 for other power-of-2 lengths, and Bluestein for arbitrary lengths.

```cpp
#include <WaveManipp.h>

std::vector<double> real = { 1.0, 0.0, -1.0, 0.0 };
std::vector<double> imag(real.size(), 0.0);

wm::dsp::fft::transform(real, imag);            // in-place forward FFT
wm::dsp::fft::inverseTransform(real, imag);     // scaled inverse FFT
```

For repeated transforms of the same length, create a reusable plan once and call it for each window. The plan caches trig tables instead of rebuilding them every transform:

```cpp
const std::size_t windowSize = 4096;
wm::dsp::fft::Plan plan(windowSize);

std::vector<double> real(windowSize);
std::vector<double> imag(windowSize);

for (/* each analysis window */) {
    // Fill real[] with windowed samples and clear imag[].
    plan.transform(real, imag);
}
```

For audio data, use `extractChannel()` when you need the complex bins, or `magnitudeSpectrum()` when a single-sided amplitude-scaled spectrum is enough:

```cpp
wm::Wave wave("audio.wav");

std::vector<double> spectrum =
    wm::dsp::fft::magnitudeSpectrum(wave, 0, 0, 4096);
```

#### [Back to Table of Contents](#table-of-contents)
