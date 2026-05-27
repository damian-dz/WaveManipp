find_package(OpenMP)

add_library(WaveManipp STATIC
    src/Wave.cpp
    src/Dsp.cpp
    src/Fft.cpp
    src/WaveBuilder.cpp
    src/WaveMixer.cpp
    src/WaveProperties.cpp
    src/RecordingBuffer.cpp
)

target_include_directories(WaveManipp PUBLIC api)

if(OpenMP_CXX_FOUND)
    target_link_libraries(WaveManipp PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
    target_compile_definitions(WaveManipp PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
