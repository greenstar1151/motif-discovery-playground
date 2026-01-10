/**
 * @file test_main.cpp
 * @brief Unit Tests for Core Operations
 * 
 * Tests each of the 5 core operations independently to verify correctness.
 */

#include "motif.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace motif;

// =============================================================================
// Test Utilities
// =============================================================================

int tests_passed = 0;
int tests_failed = 0;

void test_pass(const std::string& name) {
    std::cout << "  [PASS] " << name << std::endl;
    tests_passed++;
}

void test_fail(const std::string& name, const std::string& reason) {
    std::cout << "  [FAIL] " << name << ": " << reason << std::endl;
    tests_failed++;
}

bool approx_equal(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

// =============================================================================
// Test 1: K-mer Counting
// =============================================================================

void test_kmer_counting() {
    std::cout << "\n=== Test: K-mer Counting ===" << std::endl;
    
    // Test ZOOPS model
    {
        SequenceList seqs = {"TATATAGG", "AATATATA", "GGGGGGGG"};
        KmerCounts counts = count_kmers(seqs, 6);
        
        // "TATATA" appears in seq 0 and seq 1, so count should be 2
        if (counts["TATATA"] == 2) {
            test_pass("ZOOPS counting - TATATA appears in 2 sequences");
        } else {
            test_fail("ZOOPS counting", 
                "Expected TATATA count = 2, got " + std::to_string(counts["TATATA"]));
        }
    }
    
    // Test masked bases are skipped
    {
        SequenceList seqs = {"TATATANN", "NNATATA"};
        KmerCounts counts = count_kmers(seqs, 6);
        
        // "TATATA" at positions with N should be skipped
        // First seq: TATATA at pos 0 is valid
        // Second seq: NNATATA - NATATA invalid, ATATAN invalid, TATANN invalid
        // Actually second seq: positions are NNATA, NATAT, ATATA... none have full valid 6-mer starting with TATATA
        // Let me reconsider: seq "NNATATA" has positions 0-5=NNATAT (has N), 1-6=NATATA (has N)
        // So only first sequence contributes
        
        if (counts["TATATA"] == 1) {
            test_pass("Masked bases skipped correctly");
        } else {
            test_fail("Masked bases", 
                "Expected TATATA count = 1, got " + std::to_string(counts["TATATA"]));
        }
    }
    
    // Test total counting (non-ZOOPS)
    {
        SequenceList seqs = {"TATATATA"};  // Contains TATATA at pos 0 and pos 2
        KmerCounts counts = count_kmers_total(seqs, 6);
        
        // Two overlapping occurrences
        if (counts["TATATA"] == 2) {
            test_pass("Total counting - counts overlapping occurrences");
        } else {
            test_fail("Total counting", 
                "Expected TATATA count = 2, got " + std::to_string(counts["TATATA"]));
        }
    }
}

// =============================================================================
// Test 2: Enrichment Score
// =============================================================================

void test_enrichment_score() {
    std::cout << "\n=== Test: Enrichment Score ===" << std::endl;
    
    // Test basic enrichment calculation
    {
        KmerCounts pos_counts = {{"TATATA", 90}, {"ACGTAC", 10}};
        KmerCounts neg_counts = {{"TATATA", 5}, {"ACGTAC", 15}};
        
        EnrichmentScores scores = calculate_enrichment_scores(
            pos_counts, neg_counts, 100, 100
        );
        
        // TATATA: 90/100 - 5/100 = 0.85
        // ACGTAC: 10/100 - 15/100 = -0.05
        
        if (approx_equal(scores["TATATA"], 0.85)) {
            test_pass("Enrichment score for TATATA = 0.85");
        } else {
            test_fail("Enrichment score", 
                "Expected 0.85, got " + std::to_string(scores["TATATA"]));
        }
        
        if (approx_equal(scores["ACGTAC"], -0.05)) {
            test_pass("Depleted score for ACGTAC = -0.05");
        } else {
            test_fail("Depleted score", 
                "Expected -0.05, got " + std::to_string(scores["ACGTAC"]));
        }
    }
    
    // Test find_best_seed
    {
        EnrichmentScores scores = {{"AAA", 0.1}, {"BBB", 0.5}, {"CCC", 0.3}};
        auto [best, score] = find_best_seed(scores);
        
        if (best == "BBB" && approx_equal(score, 0.5)) {
            test_pass("find_best_seed returns highest scoring k-mer");
        } else {
            test_fail("find_best_seed", "Did not return BBB with score 0.5");
        }
    }
}

// =============================================================================
// Test 3: PWM Construction
// =============================================================================

void test_pwm_construction() {
    std::cout << "\n=== Test: PWM Construction ===" << std::endl;
    
    // Test PWM from aligned sites
    {
        std::vector<std::string> sites = {"TATATA", "TATATA", "TATATA", "TACATA"};
        PWM pwm = build_pwm_from_sites(sites, 0.0);  // No pseudocount for exact test
        
        // Position 0: all T (4/4 = 1.0)
        // Position 1: all A (4/4 = 1.0)
        // Position 2: 3T + 1C, so T should be 0.75
        
        if (pwm.width == 6) {
            test_pass("PWM width = 6");
        } else {
            test_fail("PWM width", "Expected 6, got " + std::to_string(pwm.width));
        }
        
        if (approx_equal(pwm.at(BASE_T, 0), 1.0)) {
            test_pass("Position 0: T = 1.0");
        } else {
            test_fail("Position 0", 
                "Expected T=1.0, got " + std::to_string(pwm.at(BASE_T, 0)));
        }
        
        if (approx_equal(pwm.at(BASE_T, 2), 0.75)) {
            test_pass("Position 2: T = 0.75 (3 of 4)");
        } else {
            test_fail("Position 2", 
                "Expected T=0.75, got " + std::to_string(pwm.at(BASE_T, 2)));
        }
    }
    
    // Test consensus extraction
    {
        std::vector<std::string> sites = {"TATATA", "TATATA"};
        PWM pwm = build_pwm_from_sites(sites);
        std::string consensus = get_consensus(pwm);
        
        if (consensus == "TATATA") {
            test_pass("Consensus sequence = TATATA");
        } else {
            test_fail("Consensus", "Expected TATATA, got " + consensus);
        }
    }
}

// =============================================================================
// Test 4: PWM Scoring
// =============================================================================

void test_pwm_scoring() {
    std::cout << "\n=== Test: PWM Scoring ===" << std::endl;
    
    // Build a strong PWM for TATATA
    std::vector<std::string> sites(100, "TATATA");  // 100 identical sites
    PWM pwm = build_pwm_from_sites(sites, 0.01);  // Small pseudocount
    
    // Score the matching k-mer
    {
        double score = score_kmer("TATATA", pwm);
        
        // With uniform background (0.25 each) and high PWM probs (~1.0),
        // score should be positive and significant
        if (score > 0.0) {
            test_pass("TATATA scores positive against TATATA PWM");
        } else {
            test_fail("TATATA scoring", 
                "Expected positive score, got " + std::to_string(score));
        }
    }
    
    // Score a non-matching k-mer
    {
        double bad_score = score_kmer("GGGGGG", pwm);
        double good_score = score_kmer("TATATA", pwm);
        
        if (good_score > bad_score) {
            test_pass("TATATA scores higher than GGGGGG");
        } else {
            test_fail("Scoring comparison", 
                "TATATA should score higher than GGGGGG");
        }
    }
    
    // Test sequence scanning
    {
        Sequence seq = "ACGTATATACGT";  // TATATA at position 3
        auto [pos, score] = find_best_match(seq, pwm);
        
        if (pos == 3) {
            test_pass("Best match found at position 3");
        } else {
            test_fail("Best match position", 
                "Expected 3, got " + std::to_string(pos));
        }
    }
}

// =============================================================================
// Test 5: Sequence Masking
// =============================================================================

void test_sequence_masking() {
    std::cout << "\n=== Test: Sequence Masking ===" << std::endl;
    
    // Test pattern masking
    {
        SequenceList seqs = {"ACGTATATACGT", "TATATATA"};
        SequenceList masked = mask_pattern(seqs, "TATATA");
        
        // First seq: TATATA at pos 3 -> ACGNNNNNNCGT
        // Second seq: TATATA at pos 0 -> NNNNNNTA
        
        if (masked[0] == "ACGNNNNNNCGT") {
            test_pass("First sequence masked correctly");
        } else {
            test_fail("First sequence masking", 
                "Expected ACGNNNNNNCGT, got " + masked[0]);
        }
    }
    
    // Test masked fraction calculation
    {
        SequenceList seqs = {"NNNNNNNNNN", "AAAAAAAAAA"};  // 10 masked, 10 not
        double fraction = masked_fraction(seqs);
        
        if (approx_equal(fraction, 0.5)) {
            test_pass("Masked fraction = 0.5");
        } else {
            test_fail("Masked fraction", 
                "Expected 0.5, got " + std::to_string(fraction));
        }
    }
    
    // Test in-place masking
    {
        SequenceList seqs = {"TATATAAA", "AAATATATA"};
        int count = mask_pattern_inplace(seqs, "TATATA");
        
        if (count == 2) {
            test_pass("In-place masking returned correct count");
        } else {
            test_fail("In-place masking count", 
                "Expected 2, got " + std::to_string(count));
        }
    }
}

// =============================================================================
// Integration Test: Full Pipeline
// =============================================================================

void test_full_pipeline() {
    std::cout << "\n=== Test: Full Pipeline (TATATA Discovery) ===" << std::endl;
    
    // Generate test data with known motif
    DataGenerator gen(12345);
    auto [pos_seqs, neg_seqs] = gen.generate_test_data(200, 50, "TATATA", 0.5);
    
    // Run discovery pipeline
    KmerCounts pos_counts = count_kmers(pos_seqs, 6);
    KmerCounts neg_counts = count_kmers(neg_seqs, 6);
    
    EnrichmentScores scores = calculate_enrichment_scores(
        pos_counts, neg_counts, 200, 200
    );
    
    auto [best_seed, best_score] = find_best_seed(scores);
    
    std::cout << "  Discovered motif: " << best_seed << std::endl;
    std::cout << "  Score: " << std::fixed << std::setprecision(4) << best_score << std::endl;
    
    if (best_seed == "TATATA") {
        test_pass("Correctly discovered TATATA as top motif!");
    } else {
        // Even if not exact match, check if TATATA is in top 5
        auto top = get_top_seeds(scores, 5);
        bool found = false;
        for (const auto& [kmer, _] : top) {
            if (kmer == "TATATA") {
                found = true;
                break;
            }
        }
        
        if (found) {
            test_pass("TATATA found in top 5 motifs");
        } else {
            test_fail("Pipeline test", 
                "TATATA not found in top 5. Best was: " + best_seed);
        }
    }
}

// =============================================================================
// Main
// =============================================================================

int main() {
    std::cout << "=============================================" << std::endl;
    std::cout << "  Motif Discovery Core Operations - Tests" << std::endl;
    std::cout << "=============================================" << std::endl;
    
    test_kmer_counting();
    test_enrichment_score();
    test_pwm_construction();
    test_pwm_scoring();
    test_sequence_masking();
    test_full_pipeline();
    
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  Results: " << tests_passed << " passed, " 
              << tests_failed << " failed" << std::endl;
    std::cout << "=============================================" << std::endl;
    
    return (tests_failed > 0) ? 1 : 0;
}
