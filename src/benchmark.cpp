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
            auto pwm = build_pwm_from_seed(best_seed, current_pos);
            
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
    
    // Parse command line for quick mode
    bool quick_mode = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--quick" || std::string(argv[i]) == "-q") {
            quick_mode = true;
        }
    }
    
    if (quick_mode) {
        std::cout << "\n[Quick Mode: Using smaller data sizes]\n";
        config.num_sequences = {1000, 5000};
        config.num_runs = 1;
    }
    
    auto total_start = high_resolution_clock::now();
    
    // Run all benchmarks
    benchmark_kmer_counting(config);
    benchmark_enrichment_scoring(config);
    benchmark_pwm_construction(config);
    benchmark_pwm_scoring(config);
    benchmark_sequence_masking(config);
    benchmark_full_pipeline(config);
    benchmark_scalability(config);
    estimate_memory_usage(config);
    
    auto total_end = high_resolution_clock::now();
    double total_time = duration_cast<seconds>(total_end - total_start).count();
    
    std::cout << std::endl;
    print_separator('*');
    std::cout << "  Benchmark Complete! Total time: " << total_time << " seconds" << std::endl;
    print_separator('*');
    std::cout << std::endl;
    
    return 0;
}
