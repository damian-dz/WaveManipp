#include "WaveFile.h"

#include <algorithm>
#include <limits>
#include <iostream>

namespace wm {

    void displayHeaderData(CombinedHeader ch, DATAHeader dh)
    {
        std::cout << "ChunkID: " << std::string(ch.riff.descriptor.id)
            .substr(0, sizeof(ch.riff.descriptor.id)) << std::endl;
        std::cout << "ChunkSize: " << ch.riff.descriptor.size << std::endl;
        std::cout << "Format: " << std::string(ch.riff.type)
            .substr(0, sizeof(ch.riff.type)) << std::endl;
        std::cout << "----------" << std::endl;
        std::cout << "Subchunk1ID: " << std::string(ch.wave.descriptor.id)
            .substr(0, sizeof(ch.wave.descriptor.id)) << std::endl;
        std::cout << "Subchunk1Size: " << ch.wave.descriptor.size << std::endl;
        std::cout << "AudioFormat: " << ch.wave.audioFormat << std::endl;
        std::cout << "NumChannels: " << ch.wave.numChannels << std::endl;
        std::cout << "SampleRate: " << ch.wave.sampleRate << std::endl;
        std::cout << "ByteRate: " << ch.wave.byteRate << std::endl;
        std::cout << "BlockAlign: " << ch.wave.blockAlign << std::endl;
        std::cout << "BitsPerSample: " << ch.wave.bitsPerSample << std::endl;
        std::cout << "----------" << std::endl;
        std::cout << "Subchunk2ID: " << std::string(dh.descriptor.id)
            .substr(0, sizeof(dh.descriptor.id)) << std::endl;
        std::cout << "Subchunk2Size: " << dh.descriptor.size << std::endl << std::endl;
    }

    WaveFile::WaveFile(uint8_t* data, uint32_t dataSize, uint16_t numChannels,
        uint32_t sampleRate, uint16_t bitsPerSample) :
        m_audioDataSize(dataSize),
        m_pAudioData(data)
    {
        generateHeader(dataSize, numChannels, sampleRate, bitsPerSample);
        m_pAudioData = data;
        displayHeaderData(m_combinedHeader, m_dataHeader);
    }

    WaveFile::WaveFile(uint32_t dataSize, uint16_t numChannels, uint32_t sampleRate,
        uint16_t bitsPerSample)
    {
        generateHeader(dataSize, numChannels, sampleRate, bitsPerSample);
        m_pAudioData = (uint8_t*)std::malloc(dataSize);
    }

    WaveFile::WaveFile(const std::vector<float>& samples, uint16_t numChannels,
        uint32_t sampleRate, uint16_t bitsPerSample)
    {
        int numSamples = (int)samples.size();
        if (numSamples % numChannels != 0) {
            throw "The number of samples doesn't match the number of channels";
        }
        m_audioDataSize = numSamples * bitsPerSample / 8;
        generateHeader(m_audioDataSize, numChannels, sampleRate, bitsPerSample);
        m_pAudioData = (uint8_t*)std::malloc(m_audioDataSize);

        if (bitsPerSample == 8) {
            for (int i = 0; i < numSamples; ++i) {
                m_pAudioData[i] = uint8_t(roundf(((samples[i] + 1.f)) * 127.5f));
            }
        } else if (bitsPerSample == 16) {
            int16_t* int16Data = (int16_t*)m_pAudioData;
            for (int i = 0; i < numSamples; ++i) {
                int16Data[i] = int16_t(samples[i] * 32768);
            }
        } else if (bitsPerSample == 24) {
            for (size_t i = 0; i < m_audioDataSize; i += 3) {
                uint8_t* pData = m_pAudioData;
                int val = int(samples[i / 3] * 2147483648);
                pData[i] = val >> 8;
                pData[i + 1] = val >> 16;
                pData[i + 2] = val >> 24;
            }
        } else if (bitsPerSample == 32) {
            std::memcpy(m_pAudioData, samples.data(), m_audioDataSize);
        }
    }

    WaveFile::WaveFile(const std::vector<float>& left, const std::vector<float>& right,
        uint32_t sampleRate, uint16_t bitsPerSample)
    {
        size_t samplesPerChannel = std::max(left.size(), right.size());
        int numSamples = (int)samplesPerChannel * 2;
        m_audioDataSize = numSamples * bitsPerSample / 8;
        generateHeader(m_audioDataSize, 2, sampleRate, bitsPerSample);
        m_pAudioData = (uint8_t*)std::malloc(m_audioDataSize);
        if (left.size() != right.size()) {
            if (bitsPerSample != 8) {
                memset(m_pAudioData, 0, m_audioDataSize);
            } else {
                memset(m_pAudioData, 128, m_audioDataSize);
            }
        }

        if (bitsPerSample == 8) {
            for (int i = 0; i < left.size(); ++i) {
                m_pAudioData[i << 1] = uint8_t(roundf(((left[i] + 1.f)) * 127.5f));
            }
            for (int i = 0; i < right.size(); ++i) {
                m_pAudioData[(i << 1) + 1] = uint8_t(roundf(((right[i] + 1.f)) * 127.5f));
            }
        } else if (bitsPerSample == 16) {
            int16_t* int16Data = (int16_t*)m_pAudioData;
            for (int i = 0; i < left.size(); ++i) {
                int16Data[i << 1] = int16_t(left[i] * 32768);
            }
            for (int i = 0; i < right.size(); ++i) {
                int16Data[(i << 1) + 1] = int16_t(right[i] * 32768);
            }
        } else if (bitsPerSample == 24) {
            uint8_t* pData = m_pAudioData;
            for (size_t i = 0; i < left.size(); ++i) {
                int val = int(left[i] * 2147483648);
                size_t j = (i * 3) << 1;
                pData[j] = val >> 8;
                pData[j + 1] = val >> 16;
                pData[j + 2] = val >> 24;
            }
            for (size_t i = 0; i < right.size(); ++i) {
                int val = int(right[i] * 2147483648);
                size_t j = (i * 3) << 1;
                pData[j + 3] = val >> 8;
                pData[j + 4] = val >> 16;
                pData[j + 5] = val >> 24;
            }
        } else if (bitsPerSample == 32) {
            float* floatData = (float*)m_pAudioData;
            for (int i = 0; i < left.size(); ++i) {
                floatData[i << 1] = left[i];
            }
            for (int i = 0; i < right.size(); ++i) {
                floatData[(i << 1) + 1] = right[i];
            }
        }
    }

    WaveFile::WaveFile(const std::string fileName) :
        m_pAudioData(NULL)
    {
        open(fileName);
    }

    WaveFile::WaveFile(const WaveFile& wav) :
        m_combinedHeader(wav.m_combinedHeader),
        m_dataHeader(wav.m_dataHeader),
        m_audioDataSize(wav.m_audioDataSize)
    {
        m_pAudioData = (uint8_t*)std::malloc(m_audioDataSize);
        std::memcpy(m_pAudioData, wav.m_pAudioData, m_audioDataSize);
    }


    WaveFile::~WaveFile()
    {
        if (m_pAudioData != NULL) {
            std:: free(m_pAudioData);
        }
    }

    void WaveFile::generateHeader(uint32_t dataSize, uint16_t numChannels, uint32_t sampleRate, uint16_t bitsPerSample)
    {
        if (bitsPerSample != 8 && bitsPerSample != 16 && bitsPerSample != 24 && bitsPerSample != 32) {
            throw "Bit depth not supported.";
        }

        std::memcpy(m_combinedHeader.riff.descriptor.id, "RIFF", 4);
        m_combinedHeader.riff.descriptor.size = dataSize + 36;
        std::memcpy(m_combinedHeader.riff.type, "WAVE", 4);

        std::memcpy(m_combinedHeader.wave.descriptor.id, "fmt ", 4);
        m_combinedHeader.wave.descriptor.size = 16;
        m_combinedHeader.wave.audioFormat = bitsPerSample != 32 ? 1 : 3;
        m_combinedHeader.wave.numChannels = numChannels;
        m_combinedHeader.wave.sampleRate = sampleRate;
        m_combinedHeader.wave.byteRate = sampleRate * numChannels * bitsPerSample / 8;
        m_combinedHeader.wave.blockAlign = numChannels * bitsPerSample / 8;
        m_combinedHeader.wave.bitsPerSample = bitsPerSample;

        std::memcpy(m_dataHeader.descriptor.id, "data", 4);
        m_dataHeader.descriptor.size = dataSize;
    }

    bool peekForId(const std::string id, FILE* file)
    {
        bool res = true;
        int pos = ftell(file);
        for (unsigned i = 0; i < id.size(); ++i) {
            int c = fgetc(file);
            if (id[i] != (char)c) {
                res = false;
                break;
            }
        }
        fseek(file, pos, 0);
        return res;
    }

    void findDataChunk(FILE *file)
    {
        int pos = ftell(file) + 2;
        fseek(file, pos, 0);
        std::string idData = "data";
        while (!peekForId(idData, file)) {
            pos += 2;
            fseek(file, pos, 0);
        }
    }

    void WaveFile::open(const std::string fileName)
    {
        FILE* file = fopen(fileName.c_str(), "rb");
        if (!file) {
            throw "Failed to open the file.";
        }
        fread(&m_combinedHeader, 1, sizeof(CombinedHeader), file);

        if (m_combinedHeader.wave.descriptor.size > 16) {
            int pos = ftell(file);
            fseek(file, pos + (m_combinedHeader.wave.descriptor.size - 16), 0);
        }
        if (!peekForId(std::string("data"), file)) {
            findDataChunk(file);
        }
        fread(&m_dataHeader, 1, sizeof(DATAHeader), file);
        displayHeaderData(m_combinedHeader, m_dataHeader);
        m_audioDataSize = m_dataHeader.descriptor.size;
        m_pAudioData = (uint8_t*)std::malloc(m_audioDataSize);
        fread(m_pAudioData, 1, m_audioDataSize, file);
        fclose(file);
    }

    void WaveFile::saveAs(const std::string fileName)
    {
        FILE* file = fopen(fileName.c_str(), "wb");
        if (!file) {
            throw "Failed to create the file.";
        }
        m_combinedHeader.riff.descriptor.size = m_audioDataSize + 36;
        m_combinedHeader.wave.descriptor.size = 16;
        fwrite(&m_combinedHeader, 1, sizeof(CombinedHeader), file);
        fwrite(&m_dataHeader, 1, sizeof(DATAHeader), file);
        fwrite(m_pAudioData, 1, m_audioDataSize, file);
        fclose(file);
    }

    void WaveFile::changeVolume(float volume, int channelNumber)
    {
        if (volume > 1.f || volume < 0.f) {
            throw "The value must be within the range 0-1 (inclusive)";
        }
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        float* tempData = new float[numSamples / step];
        float max = std::numeric_limits<float>::min();
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* int16Data = (int16_t*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                float res = int16Data[i] / 32768.f;
                if (abs(res) > max) {
                    max = abs(res);
                }
                tempData[i / step] = res;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 24) {
            uint8_t* pData = m_pAudioData;
            for (size_t i = offset; i < m_audioDataSize; i += step * 3) {
                float res = (pData[i] << 8 | pData[i + 1] << 16 | pData[i + 2] << 24) / 2147483648.f;
                if (abs(res) > max) {
                    max = abs(res);
                }
                tempData[i / step / 3] = res;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 32) {
            float* floatData = (float*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                float res = floatData[i];
                if (abs(res) > max) {
                    max = abs(res);
                }
                tempData[i / step] = res;
            }
        }
        float factor = volume / max;
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* int16Data = (int16_t*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                float res = tempData[i / step] * factor;
                int16Data[i] = (int16_t)(res * 32768);
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 24) {
            uint8_t* pData = m_pAudioData;
            for (size_t i = offset; i < m_audioDataSize; i += step * 3) {
                float res = tempData[i / step / 3] * factor;
                int val = int(res * 2147483648);
                pData[i] = val >> 8;
                pData[i + 1] = val >> 16;
                pData[i + 2] = val >> 24;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 32) {
            float* floatData = (float*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                floatData[i] = tempData[i / step] * factor;
            }
        }
        delete[] tempData;
    }

    std::vector<float> WaveFile::getAudioDataAsFloat(int channelNumber)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        std::vector<float> result(numSamples / step);
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* pInt16Data = (int16_t*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                result[i / step] = pInt16Data[i] / 32768.f;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 24) {
            for (size_t i = offset; i < m_audioDataSize; i += step * 3) {
                uint8_t* pData = m_pAudioData;
                result[i / step / 3] = (pData[i] << 8 | pData[i + 1] << 16 | pData[i + 2] << 24) / 2147483648.f;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 32) {
            float* floatData = (float*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                result[i / step] = floatData[i];
            }
        }
        return result;
    }

    std::vector<int32_t> WaveFile::getAudioDataAsInt(int channelNumber)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        std::vector<int32_t> result(numSamples / step);
        if (m_combinedHeader.wave.bitsPerSample == 24) {
            step *= 3;
            for (size_t i = offset * 3; i < m_audioDataSize; i += step) {
                uint8_t* pData = m_pAudioData;
                result[i / step] = (pData[i] << 8 | pData[i + 1] << 16 | pData[i + 2] << 24) >> 8;
            }
        }
        return result;
    }

    void WaveFile::setAudioFromFloatData(const std::vector<float>& data, int channelNumber)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* pInt16Data = (int16_t*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                pInt16Data[i] = (int16_t)(data[i / step] * 32768);
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 24) {
            step *= 3;
            for (size_t i = offset; i < m_audioDataSize; i += step) {
                uint8_t* pData = m_pAudioData;
                int val = int(data[i / step] * 2147483648);
                pData[i] = val >> 8;
                pData[i + 1] = val >> 16;
                pData[i + 2] = val >> 24;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 32) {
            float* pFloatData = (float*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                pFloatData[i] = data[i / step];
            }
        }
    }

    void WaveFile::setAudioFromIntData(const std::vector<int32_t>& data, int channelNumber)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        if (m_combinedHeader.wave.bitsPerSample == 24) {
            step *= 3;
            for (size_t i = offset * 3; i < m_audioDataSize; i += step) {
                uint8_t* pData = m_pAudioData;
                int32_t val = data[i / step] << 8;
                pData[i] = val >> 8;
                pData[i + 1] = val >> 16;
                pData[i + 2] = val >> 24;
            }
        }
    }

    std::vector<float> WaveFile::getAveragedOutAudioData(int binLength, int channelNumber)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        int resSize = numSamples / step / binLength;
        int remainder = (numSamples / step) % binLength;
        if (remainder != 0) {
            resSize++;
        }
        std::vector<float> result(resSize);
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* int16Data = (int16_t*)m_pAudioData;
            for (int i = 0; i < resSize; ++i) {
                double sum = 0.;
                int endPoint = 0;
                if (i == resSize - 1) {
                    endPoint = remainder == 0 ? binLength * step : remainder * step;
                } else {
                    endPoint = binLength * step;
                }
                for (int j = offset; j < endPoint; j += step) {
                    sum += int16Data[i * binLength * step + j] / 32768.f;
                }
                result[i] = (float)sum / binLength;
            }
        }
        return result;
    }

    void WaveFile::normalize(float newMin, float newMax, int channelNumber)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int offset = channelNumber;
        int step = m_combinedHeader.wave.numChannels;
        float* tempData = new float[numSamples / step];
        float min = std::numeric_limits<float>::max();
        float max = std::numeric_limits<float>::min();
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* int16Data = (int16_t*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                float res = int16Data[i] / 32768.f;
                if (res < min) {
                    min = res;
                } else if (res > max) {
                    max = res;
                }
                tempData[i / step] = res;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 24) {
            for (size_t i = offset; i < m_audioDataSize; i += step * 3) {
                uint8_t* pData = m_pAudioData;
                float res = (pData[i] << 8 | pData[i + 1] << 16 | pData[i + 2] << 24) / 2147483648.f;
                if (res < min) {
                    min = res;
                } else if (res > max) {
                    max = res;
                }
                tempData[i / step / 3] = res;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 32) {
            float* pFloatData = (float*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                float res = pFloatData[i];
                if (res < min) {
                    min = res;
                } else if (res > max) {
                    max = res;
                }
                tempData[i / step] = res;
            }
        }
        float newRange = newMax - newMin;
        float oldRange = max - min;
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* pInt16Data = (int16_t*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                float res = ((tempData[i / step] - min) * newRange) / oldRange + newMin;
                pInt16Data[i] = (int16_t)(res * 32768);
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 24) {
            for (size_t i = offset; i < m_audioDataSize; i += step * 3) {
                uint8_t* pData = m_pAudioData;
                float res = ((tempData[i / step / 3] - min) * newRange) / oldRange + newMin;
                int val = int(res * 2147483648);
                pData[i] = val >> 8;
                pData[i + 1] = val >> 16;
                pData[i + 2] = val >> 24;
            }
        } else if (m_combinedHeader.wave.bitsPerSample == 32) {
            float* pFloatData = (float*)m_pAudioData;
            for (int i = offset; i < numSamples; i += step) {
                pFloatData[i] = ((tempData[i / step] - min) * newRange) / oldRange + newMin;
            }
        }
        delete[] tempData;
    }

    uint16_t WaveFile::numChannels() const
    {
        return m_combinedHeader.wave.numChannels;
    }

    uint32_t WaveFile::sampleRate() const
    {
        return m_combinedHeader.wave.sampleRate;
    }

    void WaveFile::setSampleRate(int sampleRate)
    {
        m_combinedHeader.wave.sampleRate = sampleRate;
    }

    void WaveFile::setSampleBitDepth(int bitDepth)
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            m_combinedHeader.wave.bitsPerSample = bitDepth;
            int16_t* pInt16Data = (int16_t*)m_pAudioData;
            if (bitDepth == 8) {
                uint8_t *pData = (uint8_t*)std::malloc(numSamples);
                for (int i = 0; i < numSamples; ++i) {
                    pData[i] = uint8_t((pInt16Data[i] + 32768) >> 8);
                }
                std:: free(m_pAudioData);
                m_pAudioData = NULL;
                m_pAudioData = pData;
                m_audioDataSize /= 2;
                m_dataHeader.descriptor.size /= 2;
            } else if (bitDepth == 24) {
                m_audioDataSize = numSamples * 3;
                uint8_t* pData = (uint8_t*)std::malloc(m_audioDataSize);
                for (size_t i = 0; i < m_audioDataSize; i += 3) {
                    int val = pInt16Data[i / 3] << 16;
                    pData[i] = val >> 8;
                    pData[i + 1] = val >> 16;
                    pData[i + 2] = val >> 24;
                }
                std:: free(m_pAudioData);
                m_pAudioData = NULL;
                m_pAudioData = pData;
                m_dataHeader.descriptor.size = uint32_t(m_dataHeader.descriptor.size * 3 / 2.);
            } else if (bitDepth == 32) {
                float *pFloatData = (float*)std::malloc(numSamples * sizeof(float));
                for (int i = 0; i < numSamples; ++i) {
                    pFloatData[i] = pInt16Data[i] / 32768.f;
                }
                std:: free(m_pAudioData);
                m_pAudioData = NULL;
                m_pAudioData = (uint8_t*)pFloatData;
                m_combinedHeader.wave.audioFormat = 3;
                m_audioDataSize *= 2;
                m_dataHeader.descriptor.size *= 2;
            }
        } if (m_combinedHeader.wave.bitsPerSample == 24) {

        }
    }

    void WaveFile::append(const WaveFile& wav)
    {
        if (m_combinedHeader.wave.bitsPerSample == wav.m_combinedHeader.wave.bitsPerSample) {
            if (m_combinedHeader.wave.numChannels == wav.m_combinedHeader.wave.numChannels) {
                uint32_t combinedSize = m_audioDataSize + wav.m_audioDataSize;
                m_pAudioData = (uint8_t*)std::realloc(m_pAudioData, combinedSize);
                std::memcpy(&m_pAudioData[m_audioDataSize], wav.m_pAudioData, wav.m_audioDataSize);
                m_audioDataSize = combinedSize;
                m_combinedHeader.riff.descriptor.size = combinedSize + 36;
                m_dataHeader.descriptor.size = combinedSize;
            } else if (m_combinedHeader.wave.numChannels == 2 && wav.m_combinedHeader.wave.numChannels == 1) {
                uint32_t combinedSize = m_audioDataSize + wav.m_audioDataSize * 2;
                m_pAudioData = (uint8_t*)std::realloc(m_pAudioData, combinedSize);
                size_t frameSize = m_combinedHeader.wave.blockAlign;
                size_t bytesPerSample = m_combinedHeader.wave.bitsPerSample / 8;
                for (size_t i0 = m_audioDataSize, i1 = 0; i0 < combinedSize; i0 += frameSize, i1 += bytesPerSample) {
                    for (size_t offset = 0; offset < 2; ++offset) {
                        for (size_t j = 0; j < bytesPerSample; ++j) {
                            size_t k = offset * bytesPerSample;
                            m_pAudioData[i0 + k + j] = wav.m_pAudioData[i1 + j];
                        }
                    }
                }
                m_audioDataSize = combinedSize;
                m_combinedHeader.riff.descriptor.size = combinedSize + 36;
                m_dataHeader.descriptor.size = combinedSize;
            } else if (m_combinedHeader.wave.numChannels == 1 && wav.m_combinedHeader.wave.numChannels == 2) {
                uint32_t combinedSize = m_audioDataSize * 2 + wav.m_audioDataSize;
                uint8_t* pData = (uint8_t*)std::malloc(combinedSize);
                size_t frameSize = m_combinedHeader.wave.blockAlign * 2;
                size_t bytesPerSample = m_combinedHeader.wave.bitsPerSample / 8;
                for (size_t i0 = 0, i1 = 0; i0 < m_audioDataSize * 2; i0 += frameSize, i1 += bytesPerSample) {
                    for (size_t offset = 0; offset < 2; ++offset) {
                        for (size_t j = 0; j < bytesPerSample; ++j) {
                            size_t k = offset * bytesPerSample;
                            pData[i0 + k + j] = m_pAudioData[i1 + j];
                        }
                    }
                }
                std::memcpy(&pData[m_audioDataSize * 2], wav.m_pAudioData, wav.m_audioDataSize);
                std:: free(m_pAudioData);
                m_pAudioData = pData;
                m_audioDataSize = combinedSize;
                m_combinedHeader.wave.numChannels *= 2;
                m_combinedHeader.wave.blockAlign *= 2;
                m_combinedHeader.wave.byteRate *= 2;
                m_combinedHeader.riff.descriptor.size = combinedSize + 36;
                m_dataHeader.descriptor.size = combinedSize;
            }
        }
    }

    void WaveFile::appendAudioFromFloatData(const std::vector<float>& data)
    {

    }

    void WaveFile::downmixToMono()
    {
        int numSamples = m_audioDataSize / (m_combinedHeader.wave.bitsPerSample / 8);
        int newSize = m_audioDataSize / m_combinedHeader.wave.numChannels;
        int channelCount = m_combinedHeader.wave.numChannels;
        if (m_combinedHeader.wave.bitsPerSample == 16) {
            int16_t* newData = (int16_t*)std::malloc(newSize);
            int16_t* oldData = (int16_t*)m_pAudioData;
            for (int i = 0; i < numSamples; i += channelCount) {
                int val = 0;
                for (int c = 0; c < channelCount; ++c) {
                    val += oldData[i + c];
                }
                newData[i / channelCount] = int16_t(roundf(val / float(channelCount)));
            }
            std:: free(m_pAudioData);
            m_pAudioData = (uint8_t*)newData;
            m_combinedHeader.wave.blockAlign /= m_combinedHeader.wave.numChannels;
            m_combinedHeader.wave.byteRate /= m_combinedHeader.wave.numChannels;
            m_combinedHeader.riff.descriptor.size = newSize + 36;
            m_dataHeader.descriptor.size = newSize;
            m_audioDataSize = newSize;
            m_combinedHeader.wave.numChannels = 1;
        }
    }
}
