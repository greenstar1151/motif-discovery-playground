# Motif Discovery - Core Operations in C++

A minimal C++ implementation of fundamental motif discovery operations, designed as a reference for hardware accelerator development.

## Overview

This library implements the core computational primitives shared by motif discovery algorithms (STREME, MEME, HOMER, etc.) in plain C++17 without external dependencies or specialized optimizations.

**Key Features:**
- Pure C++ standard library implementation
- Modular, independently testable operations
- Clean abstractions suitable for hardware mapping
- Synthetic data generation for validation
- Benchmark support using both synthetic and real DNA data

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

The benchmark executable supports both synthetic and real DNA data.

### Synthetic Data (default)

```bash
# Quick sanity check
./build/motif_benchmark --quick

# Large-scale synthetic benchmark
./build/motif_benchmark --large

# Custom configuration
./build/motif_benchmark --n 50000 --l 300 --k 8 --runs 1
```

### Real DNA Data

```bash
# Run with real data (DNA1, DNA2, DNA3 datasets)
./build/motif_benchmark --real

# Run with real data for specific k-mer lengths
./build/motif_benchmark --real --k 6
./build/motif_benchmark --real --k 8

# Specify custom data files
./build/motif_benchmark --real --baseline data/upstream5000.fa --sites data/MA0007.2.sites
```

Key options:
- `--quick`: fast sanity-check with small synthetic data
- `--large`: large-scale synthetic benchmark
- `--real`: use real DNA data
- `--baseline <file>`: baseline FASTA file (default: data/upstream5000.fa)
- `--sites <file>`: binding sites file (default: data/MA0007.2.sites)
- `--n`, `--l`, `--k`: override data dimensions
- `--runs`: number of runs per configuration

## Project Structure

```
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
│       ├── data_generator.hpp    # Synthetic data generation
│       └── fasta_reader.hpp      # FASTA/.sites file I/O + injection
├── src/
│   ├── main.cpp                  # Demo program
│   ├── benchmark.cpp             # Performance benchmark
│   └── *.cpp                     # Operation implementations
├── data/                         # Real DNA data files
│   ├── upstream5000.fa           # RefSeq upstream sequences
│   └── MA0007.2.sites            # JASPAR binding sites
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

## Using Real DNA Data

This library supports loading real DNA sequences and binding sites for realistic benchmarking.

### Quick Start with Real Data

```cpp
#include "motif.hpp"
using namespace motif;

// One-liner: generate test data
auto [primary, control] = generate_real_test_data(
    "data/upstream5000.fa",   // Baseline sequences
    "data/MA0007.2.sites",    // Real binding sites
    32768,                     // Number of sequences
    1000                       // Segment length
);

// Run motif discovery pipeline
auto pos_counts = count_kmers(primary, 6);
auto neg_counts = count_kmers(control, 6);
// ...
```

### Step-by-Step Loading

```cpp
// Load baseline sequences from FASTA file
auto baseline = read_fasta("data/upstream5000.fa");
std::cout << "Loaded " << baseline.sequences.size() << " sequences\n";

// Chop 5000bp sequences into 1000bp segments
auto chopped = chop_sequences(baseline.sequences, 1000);

// Load real binding sites (JASPAR .sites format = FASTA)
auto sites = read_fasta("data/MA0007.2.sites");
std::cout << "Loaded " << sites.sequences.size() << " binding sites\n";

// Sample sequences for primary and control sets
auto primary_baseline = sample_sequences(chopped, 32768, /*seed=*/42);
auto sampled_sites = sample_sequences(sites.sequences, 32768, /*seed=*/43);

// Inject real binding sites into baseline (one per sequence)
auto primary = inject_binding_sites(primary_baseline, sampled_sites);

// Control set: pure baseline without motifs
auto control = sample_sequences(chopped, 32768, /*seed=*/1000);

// Check base composition (real DNA is not uniform!)
auto comp = compute_base_composition(primary);
std::cout << "GC content: " << (comp.gc_content() * 100) << "%\n";
```

## Morbius Benchmark Evaluation

This library includes tools for comparing motif discovery accuracy with the Morbius benchmark (IEEE e-Science 2024).

### Generating Result TSV

The `morbius_eval` tool processes FASTA datasets and outputs predictions in TSV format compatible with Morbius ground truth files.

```bash
# Basic usage
./build/morbius_eval DATASET_DNA_3.fasta DATASET_DNA_3_result.tsv

# Specify k-mer length
./build/morbius_eval DATASET_DNA_3.fasta DATASET_DNA_3_result.tsv -k 6

# Use specific seed motif
./build/morbius_eval DATASET_DNA_3.fasta DATASET_DNA_3_result.tsv --seed AGAACA

# Use control sequences instead of shuffling
./build/morbius_eval DATASET_DNA_3.fasta result.tsv --control control.fasta

# Show top N candidate seeds
./build/morbius_eval DATASET_DNA_3.fasta result.tsv --top 5
```

**Output format (TSV):**
```
seq_id	position	motif
0	859	TATATA
1	729	AGAACA
2	145	TGTACA
...
```

### Evaluating Accuracy

Compare predictions with ground truth using the Python evaluation script:

```bash
# Install dependencies
pip install pandas numpy

# Basic accuracy evaluation
python scripts/evaluate_accuracy.py DATASET_DNA_3_gt.tsv DATASET_DNA_3_result.tsv

# With position tolerance (±10bp)
python scripts/evaluate_accuracy.py DATASET_DNA_3_gt.tsv result.tsv --tolerance 10

# Save detailed per-sequence comparison
python scripts/evaluate_accuracy.py DATASET_DNA_3_gt.tsv result.tsv --output detailed.csv

# Output as JSON
python scripts/evaluate_accuracy.py DATASET_DNA_3_gt.tsv result.tsv --json
```

**Accuracy Metrics:**
- **Exact match**: Predicted position exactly matches GT position
- **Overlap match**: Predicted motif overlaps with GT motif region
- **Tolerance match**: Position within ±N bp of GT

**Example output:**
```
============================================================
  Morbius Benchmark Accuracy Evaluation
============================================================

[Dataset Statistics]
  Total sequences      : 131,072
  Predicted sequences  : 131,072
  Missing predictions  : 0
  GT motif length      : 115 bp

[Accuracy Metrics]
  Exact match          : 12,345 / 131,072
                       : 9.42%

  Overlap match        : 98,765 / 131,072
  (within GT region)   : 75.35%

[Position Error Statistics]
  Mean error           : 42.15 bp
  Median error         : 28.00 bp
  Std dev              : 35.67 bp
  Min / Max            : 0 / 884 bp

============================================================
```
