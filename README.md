# WaveManipp

**WaveManipp** is intended to be a C++ library for efficient uncompressed audio data manipulation. At this point, it supports WAV files only up to two channels. It does not support audio recording or playback on its own, and there are no plans to implement that. However, in the future, some functionalities should be incorporated into it to facilitate the cooperation with third-party playback/recording libraries.

## Introduction

In order to use the library classes and methods in our code, we have to include the necessary headers. This can be done by including the `WaveManipp.h` header file somewhere at the top of the header/source file.

```c++
#include <WaveManipp.h>
```

The easiest way to load an uncompressed audio file into memory is to do it from within the constructor by providing the path to the file (filename). All of the **WaveManipp** classes reside in the `wm` namespace.

```c++
wm::Wave wav("D:/AudioSamples/test.wav");
```

Alternatively, one can load a file by calling the empty constructor first and then using `open()`.

```c++
wm::Wave wave;
wave.open("D:/AudioSamples/test.wav");
```

The overloaded `std::ostream` operator lets us display the basic information about the opened file fairly easily.

```c++
std::cout << wave << std::endl;
```

The line above should ouput something along the lines of:

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

For example, we can now generate a sine wave and append it to the existing `Wave` object.

```c++
auto sineWave = wm::Wave::generateSine(220.f, 0.f, 44100, 44100);
wave.append(sineWave);
```

We should be able to see the updated header data.

```c++
std::cout << wave << std::endl;
```

```
ChunkID: RIFF
ChunkSize: 528248
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
Subchunk2Size: 528212
```

We can easily save the result by calling the `save()` method with a filename.

```c++
wave.save("result.wav");
```
