/**
 * @file benchmark.cpp
 * @brief Performance Benchmark for Core Motif Discovery Operations
 * 
 * Tests each operation with varying data sizes to measure performance
 * characteristics. Useful for hardware accelerator design reference.
 */

#include "motif.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include <string>

using namespace motif;
using namespace std::chrono;

// =============================================================================
// Benchmark Utilities
// =============================================================================

/**
 * @brief Timer utility for measuring execution time
 */
class Timer {
public:
    void start() { start_ = high_resolution_clock::now(); }
    void stop() { end_ = high_resolution_clock::now(); }
    
    double elapsed_ms() const {
        return duration_cast<microseconds>(end_ - start_).count() / 1000.0;
    }
    
    double elapsed_sec() const {
        return elapsed_ms() / 1000.0;
    }

private:
    time_point<high_resolution_clock> start_, end_;
};

/**
 * @brief Statistics for multiple runs
 */
struct BenchmarkStats {
    double mean_ms = 0;
    double min_ms = 0;
    double max_ms = 0;
    double std_ms = 0;
    int runs = 0;
};

BenchmarkStats compute_stats(const std::vector<double>& times) {
    BenchmarkStats stats;
    stats.runs = static_cast<int>(times.size());
    
    if (times.empty()) return stats;
    
    stats.min_ms = *std::min_element(times.begin(), times.end());
    stats.max_ms = *std::max_element(times.begin(), times.end());
    stats.mean_ms = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    
    double sq_sum = 0;
    for (double t : times) {
        sq_sum += (t - stats.mean_ms) * (t - stats.mean_ms);
    }
    stats.std_ms = std::sqrt(sq_sum / times.size());
    
    return stats;
}

void print_separator(char c = '=', int width = 80) {
    std::cout << std::string(width, c) << std::endl;
}

void print_header(const std::string& title) {
    std::cout << std::endl;
    print_separator('=');
    std::cout << "  " << title << std::endl;
    print_separator('=');
}

void print_subheader(const std::string& title) {
    std::cout << std::endl;
    std::cout << "--- " << title << " ---" << std::endl;
}

// =============================================================================
// Benchmark Configuration
// =============================================================================

struct BenchmarkConfig {
    // Data sizes to test
    std::vector<size_t> num_sequences = {1000, 5000, 10000, 20000};
    std::vector<size_t> seq_lengths = {50, 100, 200, 500};
    std::vector<size_t> kmer_lengths = {6, 8, 10};
    
    // Test parameters
    std::string target_motif = "TATATA";
    double injection_rate = 0.5;
    int num_runs = 3;  // Runs per configuration for averaging
    unsigned int seed = 42;
    
    // Real data mode
    bool use_real_data = false;
    std::string baseline_file = "data/upstream5000.fa";
    std::string sites_file = "data/MA0007.2.sites";
    size_t segment_length = 1000;
};

// =============================================================================
// Individual Operation Benchmarks
// =============================================================================

void benchmark_kmer_counting(const BenchmarkConfig& config) {
    print_header("Benchmark 1: K-mer Counting");
    
    std::cout << std::endl;
    std::cout << std::setw(12) << "Sequences" 
              << std::setw(12) << "SeqLen"
              << std::setw(8) << "K"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "K-mers Found"
              << std::setw(18) << "Throughput (seq/s)"
              << std::endl;
    print_separator('-');
    
    Timer timer;
    DataGenerator gen(config.seed);
    
    for (size_t num_seq : config.num_sequences) {
        for (size_t seq_len : {100, 200}) {  // Subset of lengths
            for (size_t k : {6, 8}) {  // Subset of k values
                // Generate data
                auto sequences = gen.generate_random_sequences(num_seq, seq_len);
                sequences = gen.inject_motif(sequences, config.target_motif, config.injection_rate);
                
                // Benchmark
                std::vector<double> times;
                KmerCounts counts;
                
                for (int run = 0; run < config.num_runs; ++run) {
                    timer.start();
                    counts = count_kmers(sequences, k);
                    timer.stop();
                    times.push_back(timer.elapsed_ms());
                }
                
                auto stats = compute_stats(times);
                double throughput = num_seq / (stats.mean_ms / 1000.0);
                
                std::cout << std::setw(12) << num_seq
                          << std::setw(12) << seq_len
                          << std::setw(8) << k
                          << std::setw(15) << std::fixed << std::setprecision(2) << stats.mean_ms
                          << std::setw(15) << counts.size()
                          << std::setw(18) << std::fixed << std::setprecision(0) << throughput
                          << std::endl;
            }
        }
    }
}

void benchmark_enrichment_scoring(const BenchmarkConfig& config) {
    print_header("Benchmark 2: Enrichment Score Calculation");
    
    std::cout << std::endl;
    std::cout << std::setw(12) << "Sequences"
              << std::setw(12) << "K-mers"
              << std::setw(15) << "Time (ms)"
              << std::setw(18) << "Best Seed"
              << std::setw(12) << "Score"
              << std::endl;
    print_separator('-');
    
    Timer timer;
    DataGenerator gen(config.seed);
    
    for (size_t num_seq : config.num_sequences) {
        size_t seq_len = 100;
        size_t k = 6;
        
        // Generate data
        auto [pos_seqs, neg_seqs] = gen.generate_test_data(
            num_seq, seq_len, config.target_motif, config.injection_rate
        );
        
        // Count k-mers first
        auto pos_counts = count_kmers(pos_seqs, k);
        auto neg_counts = count_kmers(neg_seqs, k);
        
        // Benchmark enrichment scoring
        std::vector<double> times;
        EnrichmentScores scores;
        
        for (int run = 0; run < config.num_runs; ++run) {
            timer.start();
            scores = calculate_enrichment_scores(
                pos_counts, neg_counts,
                static_cast<int>(num_seq),
                static_cast<int>(num_seq)
            );
            timer.stop();
            times.push_back(timer.elapsed_ms());
        }
        
        auto [best_seed, best_score] = find_best_seed(scores);
        auto stats = compute_stats(times);
        
        std::cout << std::setw(12) << num_seq
                  << std::setw(12) << scores.size()
                  << std::setw(15) << std::fixed << std::setprecision(3) << stats.mean_ms
                  << std::setw(18) << best_seed
                  << std::setw(12) << std::fixed << std::setprecision(4) << best_score
                  << std::endl;
    }
}

void benchmark_pwm_construction(const BenchmarkConfig& config) {
    print_header("Benchmark 3: PWM Construction");
    
    std::cout << std::endl;
    std::cout << std::setw(12) << "Sequences"
              << std::setw(12) << "Matches"
              << std::setw(15) << "Time (ms)"
              << std::setw(18) << "Consensus"
              << std::endl;
    print_separator('-');
    
    Timer timer;
    DataGenerator gen(config.seed);
    
    for (size_t num_seq : config.num_sequences) {
        // Generate data with high injection rate for many matches
        auto sequences = gen.generate_random_sequences(num_seq, 100);
        sequences = gen.inject_motif(sequences, config.target_motif, 0.8);
        
        // Benchmark PWM construction
        std::vector<double> times;
        PWM pwm;
        
        for (int run = 0; run < config.num_runs; ++run) {
            timer.start();
            pwm = build_pwm_from_seed(config.target_motif, sequences);
            timer.stop();
            times.push_back(timer.elapsed_ms());
        }
        
        // Count matches for reference
        int match_count = 0;
        for (const auto& seq : sequences) {
            if (seq.find(config.target_motif) != std::string::npos) {
                match_count++;
            }
        }
        
        auto stats = compute_stats(times);
        
        std::cout << std::setw(12) << num_seq
                  << std::setw(12) << match_count
                  << std::setw(15) << std::fixed << std::setprecision(3) << stats.mean_ms
                  << std::setw(18) << get_consensus(pwm)
                  << std::endl;
    }
}

void benchmark_pwm_scoring(const BenchmarkConfig& config) {
    print_header("Benchmark 4: PWM Scoring (Sequence Scanning)");
    
    std::cout << std::endl;
    std::cout << std::setw(12) << "Sequences"
              << std::setw(12) << "SeqLen"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Hits"
              << std::setw(18) << "Throughput (seq/s)"
              << std::endl;
    print_separator('-');
    
    Timer timer;
    DataGenerator gen(config.seed);
    
    // Build a reference PWM
    std::vector<std::string> sites(100, config.target_motif);
    PWM pwm = build_pwm_from_sites(sites, 0.1);
    
    for (size_t num_seq : config.num_sequences) {
        for (size_t seq_len : {100, 200, 500}) {
            // Generate data
            auto sequences = gen.generate_random_sequences(num_seq, seq_len);
            sequences = gen.inject_motif(sequences, config.target_motif, config.injection_rate);
            
            // Benchmark PWM scoring
            std::vector<double> times;
            int hit_count = 0;
            
            for (int run = 0; run < config.num_runs; ++run) {
                timer.start();
                hit_count = count_hits(sequences, pwm, 5.0);  // threshold = 5 bits
                timer.stop();
                times.push_back(timer.elapsed_ms());
            }
            
            auto stats = compute_stats(times);
            double throughput = num_seq / (stats.mean_ms / 1000.0);
            
            std::cout << std::setw(12) << num_seq
                      << std::setw(12) << seq_len
                      << std::setw(15) << std::fixed << std::setprecision(2) << stats.mean_ms
                      << std::setw(15) << hit_count
                      << std::setw(18) << std::fixed << std::setprecision(0) << throughput
                      << std::endl;
        }
    }
}

void benchmark_sequence_masking(const BenchmarkConfig& config) {
    print_header("Benchmark 5: Sequence Masking");
    
    std::cout << std::endl;
    std::cout << std::setw(12) << "Sequences"
              << std::setw(12) << "SeqLen"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Masked %"
              << std::setw(18) << "Throughput (seq/s)"
              << std::endl;
    print_separator('-');
    
    Timer timer;
    DataGenerator gen(config.seed);
    
    for (size_t num_seq : config.num_sequences) {
        for (size_t seq_len : {100, 200, 500}) {
            // Generate data
            auto sequences = gen.generate_random_sequences(num_seq, seq_len);
            sequences = gen.inject_motif(sequences, config.target_motif, config.injection_rate);
            
            // Benchmark masking
            std::vector<double> times;
            SequenceList masked;
            
            for (int run = 0; run < config.num_runs; ++run) {
                timer.start();
                masked = mask_pattern(sequences, config.target_motif);
                timer.stop();
                times.push_back(timer.elapsed_ms());
            }
            
            auto stats = compute_stats(times);
            double mask_pct = masked_fraction(masked) * 100.0;
            double throughput = num_seq / (stats.mean_ms / 1000.0);
            
            std::cout << std::setw(12) << num_seq
                      << std::setw(12) << seq_len
                      << std::setw(15) << std::fixed << std::setprecision(2) << stats.mean_ms
                      << std::setw(15) << std::fixed << std::setprecision(2) << mask_pct
                      << std::setw(18) << std::fixed << std::setprecision(0) << throughput
                      << std::endl;
        }
    }
}

// =============================================================================
// Full Pipeline Benchmark
// =============================================================================

void benchmark_full_pipeline(const BenchmarkConfig& config) {
    print_header("Benchmark 6: Full Discovery Pipeline (Multi-round)");
    
    std::cout << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  - Target Motif    : " << config.target_motif << std::endl;
    std::cout << "  - Injection Rate  : " << (config.injection_rate * 100) << "%" << std::endl;
    std::cout << "  - K-mer Length    : 6" << std::endl;
    std::cout << "  - Max Rounds      : 3" << std::endl;
    std::cout << std::endl;
    
    std::cout << std::setw(12) << "Sequences"
              << std::setw(12) << "SeqLen"
              << std::setw(15) << "Total (ms)"
              << std::setw(15) << "Round 1 (ms)"
              << std::setw(15) << "Found Motif"
              << std::endl;
    print_separator('-');
    
    Timer timer, round_timer;
    DataGenerator gen(config.seed);
    
    // Test with larger configurations
    std::vector<std::pair<size_t, size_t>> configs = {
        {5000, 100},
        {10000, 100},
        {10000, 200},
        {20000, 100},
        {20000, 200},
        {50000, 100},
    };
    
    for (auto [num_seq, seq_len] : configs) {
        const size_t k = 6;
        const int max_rounds = 3;
        const double threshold = 0.05;
        
        // Generate data
        auto [pos_seqs, neg_seqs] = gen.generate_test_data(
            num_seq, seq_len, config.target_motif, config.injection_rate
        );
        
        std::string discovered_motif;
        double round1_time = 0;
        
        timer.start();
        
        // Run pipeline
        SequenceList current_pos = pos_seqs;
        SequenceList current_neg = neg_seqs;
        
        for (int round = 0; round < max_rounds; ++round) {
            if (round == 0) round_timer.start();
            
            // Step 1: K-mer counting
            auto pos_counts = count_kmers(current_pos, k);
            auto neg_counts = count_kmers(current_neg, k);
            
            if (pos_counts.empty()) break;
            
            // Step 2: Enrichment scoring
            auto scores = calculate_enrichment_scores(
                pos_counts, neg_counts,
                static_cast<int>(num_seq),
                static_cast<int>(num_seq)
            );
            
            // Step 3: Find best seed
            auto [best_seed, best_score] = find_best_seed(scores);
            
            if (round == 0) {
                round_timer.stop();
                round1_time = round_timer.elapsed_ms();
                discovered_motif = best_seed;
            }
            
            if (best_score < threshold) break;
            
            // Step 4: Build PWM
            (void)build_pwm_from_seed(best_seed, current_pos);
            
            // Step 5: Mask and continue
            current_pos = mask_pattern(current_pos, best_seed);
            current_neg = mask_pattern(current_neg, best_seed);
        }
        
        timer.stop();
        
        std::cout << std::setw(12) << num_seq
                  << std::setw(12) << seq_len
                  << std::setw(15) << std::fixed << std::setprecision(1) << timer.elapsed_ms()
                  << std::setw(15) << std::fixed << std::setprecision(1) << round1_time
                  << std::setw(15) << discovered_motif
                  << std::endl;
    }
}

// =============================================================================
// Scalability Analysis
// =============================================================================

void benchmark_scalability(const BenchmarkConfig& config) {
    print_header("Benchmark 7: Scalability Analysis (K-mer Counting)");
    
    std::cout << std::endl;
    std::cout << "Testing how execution time scales with data size..." << std::endl;
    std::cout << std::endl;
    
    std::cout << std::setw(15) << "Total Bases"
              << std::setw(15) << "Time (ms)"
              << std::setw(20) << "Time/Million Bases"
              << std::endl;
    print_separator('-');
    
    Timer timer;
    DataGenerator gen(config.seed);
    
    // Vary total data size
    std::vector<std::pair<size_t, size_t>> sizes = {
        {1000, 100},    // 100K bases
        {2000, 100},    // 200K bases
        {5000, 100},    // 500K bases
        {10000, 100},   // 1M bases
        {20000, 100},   // 2M bases
        {10000, 200},   // 2M bases (different shape)
        {50000, 100},   // 5M bases
        {100000, 100},  // 10M bases
    };
    
    for (auto [num_seq, seq_len] : sizes) {
        size_t total_bases = num_seq * seq_len;
        
        auto sequences = gen.generate_random_sequences(num_seq, seq_len);
        
        std::vector<double> times;
        for (int run = 0; run < 3; ++run) {
            timer.start();
            auto counts = count_kmers(sequences, 6);
            timer.stop();
            times.push_back(timer.elapsed_ms());
        }
        
        auto stats = compute_stats(times);
        double time_per_million = stats.mean_ms / (total_bases / 1000000.0);
        
        std::cout << std::setw(15) << total_bases
                  << std::setw(15) << std::fixed << std::setprecision(1) << stats.mean_ms
                  << std::setw(20) << std::fixed << std::setprecision(2) << time_per_million
                  << std::endl;
    }
}

// =============================================================================
// Real Data Benchmark
// =============================================================================

void benchmark_real_data(const BenchmarkConfig& config) {
    print_header("Benchmark: Real DNA Data");
    
    std::cout << std::endl;
    std::cout << "Loading real data files..." << std::endl;
    std::cout << "  Baseline: " << config.baseline_file << std::endl;
    std::cout << "  Sites:    " << config.sites_file << std::endl;
    std::cout << std::endl;
    
    Timer timer, total_timer;
    total_timer.start();
    
    // Load baseline sequences
    timer.start();
    FastaData baseline_data;
    try {
        baseline_data = read_fasta(config.baseline_file);
    } catch (const std::exception& e) {
        std::cerr << "Error loading baseline file: " << e.what() << std::endl;
        std::cerr << "Make sure data files exist in data/ directory" << std::endl;
        return;
    }
    timer.stop();
    std::cout << "  Loaded " << baseline_data.sequences.size() << " baseline sequences in "
              << std::fixed << std::setprecision(1) << timer.elapsed_ms() << " ms" << std::endl;
    
    // Chop into segments
    timer.start();
    auto chopped = chop_sequences(baseline_data.sequences, config.segment_length);
    timer.stop();
    std::cout << "  Chopped into " << chopped.size() << " x " << config.segment_length 
              << "bp segments in " << std::fixed << std::setprecision(1) << timer.elapsed_ms() << " ms" << std::endl;
    
    // Load binding sites
    timer.start();
    FastaData sites_data;
    try {
        sites_data = read_fasta(config.sites_file);
    } catch (const std::exception& e) {
        std::cerr << "Error loading sites file: " << e.what() << std::endl;
        return;
    }
    timer.stop();
    std::cout << "  Loaded " << sites_data.sequences.size() << " binding sites in "
              << std::fixed << std::setprecision(1) << timer.elapsed_ms() << " ms" << std::endl;
    
    // Compute base composition
    auto comp = compute_base_composition(chopped);
    std::cout << "  GC content: " << std::fixed << std::setprecision(1) 
              << (comp.gc_content() * 100) << "%" << std::endl;
    std::cout << std::endl;
    
    // Dataset configurations (from Morbius, doi: 10.1109/e-Science62913.2024.10678700)
    std::vector<std::tuple<std::string, size_t, size_t>> datasets = {
        {"DNA1", 32768, 1024},
        {"DNA2", 65536, 2048},
        {"DNA3", 131072, 4096},
    };
    
    std::cout << std::setw(10) << "Dataset"
              << std::setw(6) << "K"
              << std::setw(12) << "Sequences"
              << std::setw(10) << "Sites"
              << std::setw(12) << "Size (MB)"
              << std::setw(15) << "K-mer (ms)"
              << std::setw(15) << "Enrich (ms)"
              << std::setw(15) << "Total (ms)"
              << std::setw(15) << "Best Seed"
              << std::endl;
    print_separator('-');
    
    for (const auto& [name, num_seq, num_sites] : datasets) {
        if (chopped.size() < num_seq) {
            std::cout << std::setw(10) << name << "  (skipped - not enough baseline sequences)" << std::endl;
            continue;
        }
        if (sites_data.sequences.size() < num_sites) {
            std::cout << std::setw(10) << name << "  (skipped - not enough binding sites)" << std::endl;
            continue;
        }
        
        // Sample baseline and sites
        auto primary_baseline = sample_sequences(chopped, num_seq, config.seed);
        auto sampled_sites = sample_sequences(sites_data.sequences, num_sites, config.seed + 1);
        
        // Inject binding sites
        timer.start();
        auto primary = inject_binding_sites(primary_baseline, sampled_sites, config.seed + 2);
        timer.stop();
        double inject_ms = timer.elapsed_ms();
        
        // Control set (different sample, no injection)
        auto control = sample_sequences(chopped, num_seq, config.seed + 1000);
        
        double size_mb = (num_seq * config.segment_length * 2) / (1024.0 * 1024.0);
        
        const std::vector<size_t> ks = config.kmer_lengths.empty()
            ? std::vector<size_t>{6}
            : config.kmer_lengths;

        for (size_t k : ks) {
            // Benchmark k-mer counting
            timer.start();
            auto pos_counts = count_kmers(primary, k);
            auto neg_counts = count_kmers(control, k);
            timer.stop();
            double kmer_ms = timer.elapsed_ms();

            // Benchmark enrichment scoring
            timer.start();
            auto scores = calculate_enrichment_scores(
                pos_counts, neg_counts,
                static_cast<int>(num_seq),
                static_cast<int>(num_seq)
            );
            auto [best_seed, best_score] = find_best_seed(scores);
            timer.stop();
            double enrich_ms = timer.elapsed_ms();

            double total_ms = inject_ms + kmer_ms + enrich_ms;

            std::cout << std::setw(10) << name
                      << std::setw(6) << k
                      << std::setw(12) << num_seq
                      << std::setw(10) << num_sites
                      << std::setw(12) << std::fixed << std::setprecision(1) << size_mb
                      << std::setw(15) << std::fixed << std::setprecision(1) << kmer_ms
                      << std::setw(15) << std::fixed << std::setprecision(2) << enrich_ms
                      << std::setw(15) << std::fixed << std::setprecision(1) << total_ms
                      << std::setw(15) << best_seed
                      << std::endl;
        }
    }
    
    total_timer.stop();
    std::cout << std::endl;
    std::cout << "Real data benchmark completed in " << std::fixed << std::setprecision(1) 
              << total_timer.elapsed_sec() << " seconds" << std::endl;
}

// =============================================================================
// Memory Estimation
// =============================================================================

void estimate_memory_usage(const BenchmarkConfig& config) {
    print_header("Memory Usage Estimation");
    
    std::cout << std::endl;
    std::cout << std::setw(12) << "Sequences"
              << std::setw(12) << "SeqLen"
              << std::setw(15) << "Seq Data (MB)"
              << std::setw(15) << "K-mers (est)"
              << std::setw(18) << "K-mer Map (MB)"
              << std::endl;
    print_separator('-');
    
    DataGenerator gen(config.seed);
    
    for (size_t num_seq : {1000, 10000, 50000, 100000}) {
        for (size_t seq_len : {100, 200}) {
            // Sequence data size
            double seq_mb = (num_seq * seq_len) / (1024.0 * 1024.0);
            
            // Estimate unique k-mers (rough approximation)
            size_t max_kmers = std::min(
                num_seq * (seq_len - 5),  // Max possible
                static_cast<size_t>(std::pow(4, 6))  // 4^6 = 4096 for k=6
            );
            
            // K-mer map size: ~(key_size + value_size + overhead) per entry
            // std::map overhead is roughly 32-48 bytes per node
            double kmer_map_mb = (max_kmers * (6 + 4 + 40)) / (1024.0 * 1024.0);
            
            std::cout << std::setw(12) << num_seq
                      << std::setw(12) << seq_len
                      << std::setw(15) << std::fixed << std::setprecision(2) << seq_mb
                      << std::setw(15) << max_kmers
                      << std::setw(18) << std::fixed << std::setprecision(2) << kmer_map_mb
                      << std::endl;
        }
    }
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    std::cout << std::endl;
    print_separator('*');
    std::cout << "  Motif Discovery Core Operations - Performance Benchmark" << std::endl;
    std::cout << "  Hardware Accelerator Design Reference" << std::endl;
    print_separator('*');
    
    BenchmarkConfig config;
    
    auto print_usage = [&]() {
        std::cout
            << "\nUsage: motif_benchmark [options]\n\n"
            << "Options:\n"
            << "  -q, --quick         Use small data sizes for fast checks\n"
            << "  --large             Use large data sizes (seconds to tens of seconds)\n"
            << "  --real              Use real DNA data\n"
            << "  --baseline <file>   Baseline FASTA file (default: data/upstream5000.fa)\n"
            << "  --sites <file>      Binding sites file (default: data/MA0007.2.sites)\n"
            << "  --n <int>           Override sequences per set\n"
            << "  --l <int>            Override sequence length\n"
            << "  --k <int>            Override k-mer length\n"
            << "  --runs <int>         Number of runs per configuration\n"
            << "  --motif <string>     Motif to inject (default TATATA)\n"
            << "  --rate <float>       Injection rate in positive set (0..1)\n"
            << "  --seed <int>         RNG seed\n"
            << "  -h, --help           Show this help\n";
    };
    
    // Parse command line
    bool quick_mode = false;
    bool large_mode = false;
    bool override_n = false;
    bool override_l = false;
    bool override_k = false;
    bool override_runs = false;
    
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                print_usage();
                std::exit(2);
            }
            return argv[++i];
        };
        
        auto parse_size_t = [&](const std::string& s, const char* name) -> size_t {
            try {
                size_t pos = 0;
                const auto v = std::stoull(s, &pos, 10);
                if (pos != s.size()) throw std::invalid_argument("bad");
                return static_cast<size_t>(v);
            } catch (...) {
                std::cerr << "Invalid value for " << name << "\n";
                print_usage();
                std::exit(2);
            }
        };
        
        auto parse_int = [&](const std::string& s, const char* name) -> int {
            try {
                size_t pos = 0;
                const int v = std::stoi(s, &pos, 10);
                if (pos != s.size()) throw std::invalid_argument("bad");
                return v;
            } catch (...) {
                std::cerr << "Invalid value for " << name << "\n";
                print_usage();
                std::exit(2);
            }
        };
        
        auto parse_double = [&](const std::string& s, const char* name) -> double {
            try {
                size_t pos = 0;
                const double v = std::stod(s, &pos);
                if (pos != s.size()) throw std::invalid_argument("bad");
                return v;
            } catch (...) {
                std::cerr << "Invalid value for " << name << "\n";
                print_usage();
                std::exit(2);
            }
        };
        
        if (arg == "--quick" || arg == "-q") {
            quick_mode = true;
        } else if (arg == "--large") {
            large_mode = true;
        } else if (arg == "--real") {
            config.use_real_data = true;
        } else if (arg == "--baseline") {
            config.baseline_file = require_value("--baseline");
        } else if (arg == "--sites") {
            config.sites_file = require_value("--sites");
        } else if (arg == "--n") {
            const auto v = parse_size_t(require_value("--n"), "--n");
            config.num_sequences = {v};
            override_n = true;
        } else if (arg == "--l") {
            const auto v = parse_size_t(require_value("--l"), "--l");
            config.seq_lengths = {v};
            config.segment_length = v;
            override_l = true;
        } else if (arg == "--k") {
            const auto v = parse_size_t(require_value("--k"), "--k");
            config.kmer_lengths = {v};
            override_k = true;
        } else if (arg == "--runs") {
            const auto v = parse_int(require_value("--runs"), "--runs");
            if (v <= 0) {
                std::cerr << "--runs must be > 0\n";
                std::exit(2);
            }
            config.num_runs = v;
            override_runs = true;
        } else if (arg == "--motif") {
            config.target_motif = require_value("--motif");
        } else if (arg == "--rate") {
            const double v = parse_double(require_value("--rate"), "--rate");
            if (v < 0.0 || v > 1.0) {
                std::cerr << "--rate must be between 0 and 1\n";
                std::exit(2);
            }
            config.injection_rate = v;
        } else if (arg == "--seed") {
            const int v = parse_int(require_value("--seed"), "--seed");
            if (v < 0) {
                std::cerr << "--seed must be >= 0\n";
                std::exit(2);
            }
            config.seed = static_cast<unsigned int>(v);
        } else if (arg == "-h" || arg == "--help") {
            print_usage();
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage();
            return 2;
        }
    }
    
    if (large_mode) {
        std::cout << "\n[Large Mode: Using larger data sizes]\n";
        config.num_sequences = {20000, 50000, 100000};
        config.seq_lengths = {200, 500};
        config.kmer_lengths = {6, 8};
        config.num_runs = 1;
    }
    
    if (quick_mode) {
        std::cout << "\n[Quick Mode: Using smaller data sizes]\n";
        config.num_sequences = {1000, 5000};
        config.seq_lengths = {50, 100};
        config.kmer_lengths = {6};
        config.num_runs = 1;
    }
    
    // If user overrides any size, keep other defaults unless explicitly overridden.
    if ((override_n || override_l || override_k) && !override_runs) {
        config.num_runs = 1;
    }
    
    auto total_start = high_resolution_clock::now();
    
    // Run benchmarks based on mode
    if (config.use_real_data) {
        // Real data mode - only run real data benchmark
        benchmark_real_data(config);
    } else {
        // Synthetic data mode - run all synthetic benchmarks
        benchmark_kmer_counting(config);
        benchmark_enrichment_scoring(config);
        benchmark_pwm_construction(config);
        benchmark_pwm_scoring(config);
        benchmark_sequence_masking(config);
        benchmark_full_pipeline(config);
        benchmark_scalability(config);
        estimate_memory_usage(config);
    }
    
    auto total_end = high_resolution_clock::now();
    double total_time = duration_cast<seconds>(total_end - total_start).count();
    
    std::cout << std::endl;
    print_separator('*');
    std::cout << "  Benchmark Complete! Total time: " << total_time << " seconds" << std::endl;
    print_separator('*');
    std::cout << std::endl;
    
    return 0;
}
