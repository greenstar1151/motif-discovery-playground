# Motif Discovery - Core Operations in C++

A minimal C++ implementation of fundamental motif discovery operations, designed as a reference for hardware accelerator development.

## Overview

This library implements the core computational primitives shared by motif discovery algorithms (STREME, MEME, HOMER, etc.) in plain C++17 without external dependencies or specialized optimizations.

**Key Features:**
- Pure C++ standard library implementation
- Modular, independently testable operations
- Clean abstractions suitable for hardware mapping
- Synthetic data generation for validation

## Core Operations

| # | Operation | File | Complexity |
|---|-----------|------|------------|
| 1 | K-mer Counting | `kmer_counting.cpp` | O(N × L) |
| 2 | Enrichment Scoring | `enrichment_score.cpp` | O(4^k) |
| 3 | PWM Construction | `pwm_construction.cpp` | O(M × W) |
| 4 | PWM Scoring | `pwm_scoring.cpp` | O(W) |
| 5 | Sequence Masking | `sequence_masking.cpp` | O(N × L) |

Where: N = sequences, L = avg length, W = motif width, M = aligned sites

## Building

```bash
cmake -B build && cmake --build build
```

## Running

```bash
# Run the demo
./build/motif_discovery

# Run unit tests
./build/motif_test

# Run benchmark (large-scale preset)
./build/motif_benchmark --large
```

## Benchmarking

The benchmark executable generates large synthetic datasets and times each
core operation (k-mer counting, enrichment, PWM build/scan, masking). Use
`--large` for a seconds-to-tens-of-seconds run, or override sizes with
`--n`, `--l`, `--k`.

```bash
# Few seconds
./build/motif_benchmark --n 20000 --l 200 --k 8 --runs 1

# Tens of seconds (heavier)
./build/motif_benchmark --large

# Single-size focused run
./build/motif_benchmark --n 50000 --l 300 --k 8 --runs 1
```

Key options:
- `--n`: sequences per set (positive/negative)
- `--l`: sequence length
- `--k`: k-mer length
- `--runs`: number of runs per configuration
- `--large`: large-scale preset
- `--quick`: fast sanity-check preset

## Project Structure

```
cpp/
├── CMakeLists.txt
├── include/
│   ├── motif.hpp                 # Main header
│   └── motif/
│       ├── types.hpp             # Common types (PWM, Sequence, etc.)
│       ├── kmer_counting.hpp
│       ├── enrichment_score.hpp
│       ├── pwm_construction.hpp
│       ├── pwm_scoring.hpp
│       ├── sequence_masking.hpp
│       └── data_generator.hpp
├── src/
│   ├── main.cpp                  # Demo program
│   └── *.cpp                     # Operation implementations
└── tests/
    └── test_main.cpp
```

## Usage Example

```cpp
#include "motif.hpp"
using namespace motif;

// Generate test data with embedded motif
DataGenerator gen(42);
auto [primary, control] = gen.generate_test_data(200, 50, "TATATA", 0.5);

// Count k-mers (ZOOPS model)
KmerCounts pos_counts = count_kmers(primary, 6);
KmerCounts neg_counts = count_kmers(control, 6);

// Calculate differential enrichment
EnrichmentScores scores = calculate_enrichment_scores(
    pos_counts, neg_counts, 200, 200
);

// Find best seed
auto [best_seed, best_score] = find_best_seed(scores);

// Build PWM from seed matches
PWM pwm = build_pwm_from_seed(best_seed, primary);

// Mask discovered motif for next iteration
primary = mask_pattern(primary, best_seed);
```

## Data Structures

### PWM (Position Weight Matrix)
```cpp
struct PWM {
    size_t width;
    double matrix[4][MAX_WIDTH];  // [A,C,G,T][position]
};
```

### K-mer Counts
```cpp
using KmerCounts = std::map<std::string, int>;
```

## References

- Bailey, T.L. (2021). "STREME: accurate and versatile sequence motif discovery." *Bioinformatics*, 37(18), 2834-2840. [DOI](https://doi.org/10.1093/bioinformatics/btab203)
- [MEME Suite Documentation](https://meme-suite.org/meme/doc/streme.html)
