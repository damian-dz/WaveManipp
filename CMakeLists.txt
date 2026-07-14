find_package(OpenMP)

add_library(WaveManipp STATIC
    src/Wave.cpp
    src/Dsp.cpp
    src/PhaseVocoder.cpp
    src/Fft.cpp
    src/WaveBuilder.cpp
    src/WaveMixer.cpp
    src/WaveProperties.cpp
    src/RecordingBuffer.cpp
    src/WaveFileFormat.cpp
    src/WaveFileReader.cpp
    src/WaveFileWriter.cpp
)

target_include_directories(WaveManipp PUBLIC api)

if(OpenMP_CXX_FOUND)
    target_link_libraries(WaveManipp PUBLIC OpenMP::OpenMP_CXX)
endif()

if(MSVC)
    target_compile_definitions(WaveManipp PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# WaveManipp's own tests: mirrors app/CMakeLists.txt's WaveformPeakCache_tests target,
# but for tests that exercise WaveManipp directly (no app/src involvement at all), so
# they don't need to link wxWidgets or pull in any app-level sources.
include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(doctest)

add_executable(WaveManipp_tests
    tests/WaveManippTestMain.cpp
    tests/PhaseVocoderTimeStretcherTests.cpp
    tests/ReverbProcessorTests.cpp
    tests/WaveRevisionTests.cpp
    tests/WaveFileReaderTests.cpp
    tests/WaveFileWriterTests.cpp
)
target_include_directories(WaveManipp_tests PRIVATE src)
target_link_libraries(WaveManipp_tests PRIVATE WaveManipp doctest::doctest)
if(MSVC)
    target_compile_definitions(WaveManipp_tests PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

include(CTest)
add_test(NAME WaveManipp_tests COMMAND WaveManipp_tests)

# Manual dev tool, not run by ctest and excluded from normal builds (EXCLUDE_FROM_ALL):
# lets an external process kill a multi-block WaveFileWriter write mid-stream, then
# verifies the last checkpoint left a valid, readable file (docs/AUDIO_STREAMING_PLAN.md's
# recoverable-prefix contract). Build explicitly with
# `cmake --build . --target WaveFileWriterInterruptDemo`.
add_executable(WaveFileWriterInterruptDemo EXCLUDE_FROM_ALL tests/manual/WaveFileWriterInterruptDemo.cpp)
target_link_libraries(WaveFileWriterInterruptDemo PRIVATE WaveManipp)
if(MSVC)
    target_compile_definitions(WaveFileWriterInterruptDemo PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
