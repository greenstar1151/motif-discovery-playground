cmake_minimum_required(VERSION 3.14)
project(MotifDiscovery 
    VERSION 1.0.0 
    LANGUAGES CXX
    DESCRIPTION "Core operations for motif discovery - hardware reference implementation"
)

# Modern C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings for clean code
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Option for building tests
option(BUILD_TESTS "Build test programs" ON)

# -----------------------------------------------------------------------------
# Core Library
# -----------------------------------------------------------------------------
add_library(motif_core STATIC
    src/kmer_counting.cpp
    src/enrichment_score.cpp
    src/pwm_construction.cpp
    src/pwm_scoring.cpp
    src/sequence_masking.cpp
    src/data_generator.cpp
)

target_include_directories(motif_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# -----------------------------------------------------------------------------
# Main Executable (Motif Discovery Demo)
# -----------------------------------------------------------------------------
add_executable(motif_discovery
    src/main.cpp
)

target_link_libraries(motif_discovery PRIVATE motif_core)

# -----------------------------------------------------------------------------
# Test Executable
# -----------------------------------------------------------------------------
if(BUILD_TESTS)
    add_executable(motif_test
        tests/test_main.cpp
    )
    target_link_libraries(motif_test PRIVATE motif_core)
endif()

# -----------------------------------------------------------------------------
# Benchmark Executable
# -----------------------------------------------------------------------------
add_executable(motif_benchmark
    src/benchmark.cpp
)

target_link_libraries(motif_benchmark PRIVATE motif_core)

# -----------------------------------------------------------------------------
# Morbius Evaluation Executable
# -----------------------------------------------------------------------------
add_executable(morbius_eval
    src/morbius_eval.cpp
)

target_link_libraries(morbius_eval PRIVATE motif_core)

# -----------------------------------------------------------------------------
# Aggressive Optimizations for Release Build
# -----------------------------------------------------------------------------
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        # Core library optimizations
        target_compile_options(motif_core PRIVATE
            -O3
            -march=native
            -mtune=native
            -funroll-loops
            -fomit-frame-pointer
            -ffast-math
            -flto
        )
        target_compile_definitions(motif_core PRIVATE NDEBUG)
        
        # Benchmark optimizations
        target_compile_options(motif_benchmark PRIVATE
            -O3
            -march=native
            -mtune=native
            -funroll-loops
            -fomit-frame-pointer
            -ffast-math
            -flto
        )
        target_link_options(motif_benchmark PRIVATE -flto)
        target_compile_definitions(motif_benchmark PRIVATE NDEBUG)
        
        # Morbius eval optimizations
        target_compile_options(morbius_eval PRIVATE
            -O3
            -march=native
            -mtune=native
            -funroll-loops
            -fomit-frame-pointer
            -ffast-math
            -flto
        )
        target_link_options(morbius_eval PRIVATE -flto)
        target_compile_definitions(morbius_eval PRIVATE NDEBUG)
    endif()
endif()

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------
install(TARGETS motif_discovery DESTINATION bin)
install(TARGETS motif_core DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Motif Discovery Configuration ===")
message(STATUS "  C++ Standard    : C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Build Tests     : ${BUILD_TESTS}")
message(STATUS "")
