/**
 * @file main.cpp
 * @brief Motif Discovery Demo - Full Pipeline
 * 
 * Demonstrates the complete motif discovery workflow using all 5 core operations.
 * Matches the behavior of motif_discovery_tutorial.py
 */

#include "motif.hpp"
#include <iostream>
#include <iomanip>

using namespace motif;

/**
 * @brief Print a separator line
 */
void print_separator(char c = '=', int width = 60) {
    std::cout << std::string(width, c) << std::endl;
}

/**
 * @brief Run the simple motif discovery simulation
 */
void run_motif_discovery() {
    // =========================================================================
    // Configuration
    // =========================================================================
    const std::string TARGET_MOTIF = "TATATA";
    const size_t NUM_SEQUENCES = 200;
    const size_t SEQ_LENGTH = 50;
    const size_t K = 6;  // k-mer length
    const double SIGNIFICANCE_THRESHOLD = 0.1;
    const int MAX_ROUNDS = 3;
    
    print_separator();
    std::cout << "       Motif Discovery - C++ Implementation" << std::endl;
    print_separator();
    std::cout << std::endl;
    
    // =========================================================================
    // Step 0: Data Preparation
    // =========================================================================
    std::cout << "[Setup]" << std::endl;
    std::cout << "  Target Hidden Motif : " << TARGET_MOTIF << std::endl;
    std::cout << "  Sequence Length     : " << SEQ_LENGTH << std::endl;
    std::cout << "  Number of Sequences : " << NUM_SEQUENCES << std::endl;
    std::cout << "  k-mer Length        : " << K << std::endl;
    std::cout << std::endl;
    
    // Generate test data
    DataGenerator generator(42);  // Fixed seed for reproducibility
    auto [pos_seqs, neg_seqs] = generator.generate_test_data(
        NUM_SEQUENCES, SEQ_LENGTH, TARGET_MOTIF, 0.5
    );
    
    std::cout << "[Data]" << std::endl;
    std::cout << "  Positive Set Size : " << pos_seqs.size() << std::endl;
    std::cout << "  Negative Set Size : " << neg_seqs.size() << std::endl;
    std::cout << std::endl;
    
    // =========================================================================
    // Iterative Discovery Loop
    // =========================================================================
    SequenceList current_pos_seqs = pos_seqs;
    SequenceList current_neg_seqs = neg_seqs;
    
    for (int round = 1; round <= MAX_ROUNDS; ++round) {
        print_separator('-');
        std::cout << "  Round " << round << std::endl;
        print_separator('-');
        
        // ---------------------------------------------------------------------
        // Step 1: K-mer Counting
        // ---------------------------------------------------------------------
        KmerCounts pos_counts = count_kmers(current_pos_seqs, K);
        KmerCounts neg_counts = count_kmers(current_neg_seqs, K);
        
        if (pos_counts.empty()) {
            std::cout << "  No k-mers found. Stopping." << std::endl;
            break;
        }
        
        // ---------------------------------------------------------------------
        // Step 2: Enrichment Scoring
        // ---------------------------------------------------------------------
        EnrichmentScores scores = calculate_enrichment_scores(
            pos_counts, neg_counts,
            static_cast<int>(pos_seqs.size()),
            static_cast<int>(neg_seqs.size())
        );
        
        // Find best seed
        auto [best_seed, best_score] = find_best_seed(scores);
        
        std::cout << "  Best Seed Found  : " << best_seed << std::endl;
        std::cout << "  Enrichment Score : " << std::fixed << std::setprecision(4) 
                  << best_score << std::endl;
        
        // Get counts for display
        int pos_count = pos_counts.count(best_seed) ? pos_counts.at(best_seed) : 0;
        int neg_count = neg_counts.count(best_seed) ? neg_counts.at(best_seed) : 0;
        std::cout << "  Count (Pos/Neg)  : " << pos_count << " / " << neg_count << std::endl;
        
        // Check significance threshold
        if (best_score < SIGNIFICANCE_THRESHOLD) {
            std::cout << std::endl;
            std::cout << "  -> Score below threshold (" << SIGNIFICANCE_THRESHOLD 
                      << "). Stopping." << std::endl;
            break;
        }
        
        // ---------------------------------------------------------------------
        // Step 3: PWM Construction
        // ---------------------------------------------------------------------
        PWM pwm = build_pwm_from_seed(best_seed, current_pos_seqs);
        std::cout << "  PWM:" << std::endl;
        print_pwm(pwm, 2);
        
        // Show consensus
        std::cout << "  Consensus: " << get_consensus(pwm) << std::endl;
        
        // ---------------------------------------------------------------------
        // Step 4: PWM Scoring (optional demonstration)
        // ---------------------------------------------------------------------
        double seed_score = score_kmer(best_seed, pwm);
        std::cout << "  Seed PWM Score: " << std::fixed << std::setprecision(2) 
                  << seed_score << " bits" << std::endl;
        
        // ---------------------------------------------------------------------
        // Step 5: Sequence Masking
        // ---------------------------------------------------------------------
        std::cout << "  Erasing '" << best_seed << "' from sequences..." << std::endl;
        current_pos_seqs = mask_pattern(current_pos_seqs, best_seed);
        current_neg_seqs = mask_pattern(current_neg_seqs, best_seed);
        
        std::cout << "  Masked fraction (pos): " << std::fixed << std::setprecision(2)
                  << (masked_fraction(current_pos_seqs) * 100) << "%" << std::endl;
        std::cout << std::endl;
    }
    
    print_separator();
    std::cout << "  Motif Discovery Complete!" << std::endl;
    print_separator();
}

int main() {
    run_motif_discovery();
    return 0;
}
