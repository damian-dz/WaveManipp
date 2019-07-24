# WaveManipp

**WaveManipp** is intended to be a library for efficient uncompressed audio data manipulation. At this point, it supports WAV files only up to two channels.

## Introduction

In order to use the library classes and methods in our code, we have to include the necessary headers. This can be done by including the `WaveManipp.h` header file somewhere at the top of the page.

```c++
#include <WaveManipp.h>
```

The easiest way to load an uncompressed audio file into memory is to do it from within the constructor by proving the path to the file (filename). All of the classes reside in the `wm` namespace.

```c++
wm::Wave wav("D:/AudioSamples/test.wav");
```

Alternatively, one can load a file by calling the empty constructor first and then using `open()`.

```c++
wm::Wave wav;
wav.open("D:/AudioSamples/test.wav");
```

The overloaded `std::ostream` operator lets us display the basic information about the opened file fairly easily.

```c++
std::cout << wav << std::endl;
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
