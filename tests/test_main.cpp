/**
 * @file test_main.cpp
 * @brief Unit Tests for Core Operations
 * 
 * Tests each of the 5 core operations independently to verify correctness.
 */

#include "motif.hpp"
#include <cassert>
#include <cmath>
#include <fstream>
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
// Test 6: Markov Background Model
// =============================================================================

void test_markov_background() {
    std::cout << "\n=== Test: Markov Background Model ===" << std::endl;

    SequenceList control(20, "AAAAAA");
    control.push_back("AAACAA");

    MarkovOrder3 bg;
    bg.train(control, 0.1);

    const double p_a_given_aaa = bg.conditional_prob("AAAAAA", 3); // P(A|AAA)
    const double p_c_given_aaa = bg.conditional_prob("AAACAA", 3); // P(C|AAA)

    if (p_a_given_aaa > p_c_given_aaa) {
        test_pass("P(A|AAA) > P(C|AAA) after training");
    } else {
        test_fail("Markov conditional probability ordering",
                  "Expected P(A|AAA) > P(C|AAA)");
    }

    const double p_out_of_range = bg.conditional_prob("AAAAAA", 999);
    if (approx_equal(p_out_of_range, 0.25)) {
        test_pass("Out-of-range conditional_prob returns 0.25");
    } else {
        test_fail("Out-of-range conditional_prob",
                  "Expected 0.25, got " + std::to_string(p_out_of_range));
    }

    const double p_invalid_base = bg.conditional_prob("AAANAA", 3);
    if (approx_equal(p_invalid_base, 0.25)) {
        test_pass("Invalid-base conditional_prob returns 0.25");
    } else {
        test_fail("Invalid-base conditional_prob",
                  "Expected 0.25, got " + std::to_string(p_invalid_base));
    }
}

// =============================================================================
// Test 7: Seed Discovery Utilities
// =============================================================================

void test_seed_discovery_utils() {
    std::cout << "\n=== Test: Seed Discovery Utilities ===" << std::endl;

    {
        const double z = z_score_two_proportion(80, 100, 10, 100);
        if (z > 0.0) {
            test_pass("Two-proportion z-score is positive for enriched primary");
        } else {
            test_fail("z_score_two_proportion", "Expected positive z-score");
        }
    }

    {
        const size_t hd = hamming_distance("AAAA", "AAAT");
        if (hd == 1) {
            test_pass("Hamming distance exact mismatch count");
        } else {
            test_fail("hamming_distance", "Expected 1, got " + std::to_string(hd));
        }
    }

    {
        const size_t hd_mismatch_len = hamming_distance("AAA", "AA");
        if (hd_mismatch_len == std::numeric_limits<size_t>::max()) {
            test_pass("Hamming distance returns SIZE_MAX for different lengths");
        } else {
            test_fail("hamming_distance length mismatch",
                      "Expected SIZE_MAX for different lengths");
        }
    }

    {
        SequenceList seqs = {
            "GGGTATATACCC",
            "TATATAAAAAAA",
            "CCCCCCCCCCCC"
        };

        const auto sites = collect_hd_sites("TATATA", seqs, 6, 1);
        if (sites.size() == 2) {
            test_pass("collect_hd_sites returns expected number of matches");
        } else {
            test_fail("collect_hd_sites count",
                      "Expected 2, got " + std::to_string(sites.size()));
        }

        const auto wide_sites = collect_hd_sites("TATATA", {"CCCTATATAGGG"}, 8, 1);
        if (wide_sites.size() == 1 && wide_sites[0].size() == 8) {
            test_pass("collect_hd_sites target_width expansion works");
        } else {
            test_fail("collect_hd_sites target_width",
                      "Expected exactly one width-8 site");
        }
    }

    {
        SequenceList primary(40, "GGGTATATACCC");
        SequenceList control(40, "GGGCGCGCGCCC");

        const auto pool = build_seed_pool(primary, control, {6}, 1.0, 10);
        const bool has_seed = std::any_of(
            pool.begin(), pool.end(),
            [](const SeedCandidate& c) { return c.kmer == "TATATA"; }
        );

        if (!pool.empty() && has_seed) {
            test_pass("build_seed_pool includes enriched motif candidate");
        } else {
            test_fail("build_seed_pool",
                      "Expected non-empty pool containing TATATA");
        }
    }
}

// =============================================================================
// Test 8: EM Refinement + Mask Utilities
// =============================================================================

void test_em_refinement_and_masking() {
    std::cout << "\n=== Test: EM Refinement + Mask Utilities ===" << std::endl;

    {
        SequenceList seqs = {"ACGTACGT"};
        MaskMatrix mask = init_mask(seqs);
        mark_window(mask[0], 2, 3);

        const SequenceList masked = apply_mask(seqs, mask);
        const double frac = masked_fraction(mask);

        if (masked[0] == "ACNNNCGT") {
            test_pass("apply_mask respects marked window");
        } else {
            test_fail("apply_mask", "Expected ACNNNCGT, got " + masked[0]);
        }

        if (window_is_masked(mask[0], 2, 3) && approx_equal(frac, 3.0 / 8.0)) {
            test_pass("window_is_masked and masked_fraction are consistent");
        } else {
            test_fail("Mask utility consistency", "Unexpected mask window/fraction behavior");
        }
    }

    DataGenerator gen(7);
    auto [primary, control] = gen.generate_test_data(300, 60, "TATATA", 0.6);

    MarkovOrder3 bg;
    bg.train(control);

    const auto init_sites = collect_hd_sites("TATATA", primary, 6, 1);
    if (init_sites.size() < 8) {
        test_fail("EM init sites", "Too few initialization sites for EM");
        return;
    }

    const PWM init_pwm = build_pwm_from_sites(init_sites, 0.2);
    EMConfig cfg;
    cfg.max_iterations = 8;
    cfg.patience = 2;
    cfg.min_improvement = 1e-3;
    cfg.min_sites = 8;
    cfg.max_mstep_sites = 20000;
    cfg.pseudocount = 0.2;

    const auto refined = run_em_refinement(init_pwm, primary, control, bg, cfg, "TATATA");
    if (refined.has_value() && refined->pwm.width == 6 && refined->enrichment > 1.0) {
        test_pass("run_em_refinement (no mask) produces enriched model");
    } else {
        test_fail("run_em_refinement (no mask)", "Expected valid enriched refined model");
    }

    const PWM motif_pwm = build_pwm_from_sites(std::vector<std::string>(30, "TATATA"), 0.1);
    const MatchResult best = find_best_match_llr("GGGTATATAGGG", motif_pwm, bg);
    if (best.position >= 0 && std::isfinite(best.score)) {
        test_pass("find_best_match_llr returns finite best match");
    } else {
        test_fail("find_best_match_llr", "Expected finite match with valid position");
    }

    MaskRow full_mask(12, 1);
    const MatchResult masked_best = find_best_match_llr("GGGTATATAGGG", full_mask, motif_pwm, bg);
    if (masked_best.position == -1) {
        test_pass("find_best_match_llr (masked) skips fully masked sequence");
    } else {
        test_fail("find_best_match_llr (masked)", "Expected no valid position under full mask");
    }

    MaskMatrix pmask = init_mask(primary);
    MaskMatrix cmask = init_mask(control);
    const auto refined_masked = run_em_refinement(
        init_pwm, primary, control, pmask, cmask, bg, cfg, "TATATA"
    );

    if (refined_masked.has_value() && refined_masked->enrichment > 1.0) {
        test_pass("run_em_refinement (with mask) produces enriched model");
    } else {
        test_fail("run_em_refinement (with mask)", "Expected valid enriched refined model");
    }

    if (refined.has_value()) {
        const int erased = erase_by_pwm(
            primary, pmask, refined->pwm, bg, refined->learned_threshold
        );
        if (erased > 0 && masked_fraction(pmask) > 0.0) {
            test_pass("erase_by_pwm masks at least one hit window");
        } else {
            test_fail("erase_by_pwm", "Expected erased > 0 and non-zero masked fraction");
        }
    }
}

// =============================================================================
// Test 9: IO Utilities
// =============================================================================

void test_io_utils() {
    std::cout << "\n=== Test: IO Utilities ===" << std::endl;

    {
        const std::string seq = "ACGTACGTACGT";
        std::mt19937 rng(123);
        const std::string shuffled = shuffle_sequence_kmer_preserving(seq, rng, 2);

        const KmerCounts original_di = count_kmers_total({seq}, 2);
        const KmerCounts shuffled_di = count_kmers_total({shuffled}, 2);

        if (shuffled.size() == seq.size()) {
            test_pass("k-mer preserving shuffle keeps sequence length");
        } else {
            test_fail("shuffle length", "Expected same sequence length after shuffle");
        }

        if (original_di == shuffled_di) {
            test_pass("Dinucleotide counts preserved after k=2 shuffle");
        } else {
            test_fail("k-mer preserving shuffle", "Dinucleotide composition changed");
        }
    }

    {
        SequenceList input = {"ACGTACGT", "TATATATA"};
        const auto control_a = generate_control_sequences(input, 42, 2);
        const auto control_b = generate_control_sequences(input, 42, 2);

        if (control_a == control_b) {
            test_pass("generate_control_sequences is deterministic for fixed seed");
        } else {
            test_fail("generate_control_sequences determinism",
                      "Expected same output with same seed");
        }
    }

    {
        std::string jaspar_path = "data/MA0007.2.jaspar";
        std::ifstream f(jaspar_path);
        if (!f.is_open()) {
            jaspar_path = "../data/MA0007.2.jaspar";
            f.clear();
            f.open(jaspar_path);
        }

        if (!f.is_open()) {
            test_fail("load_jaspar_pwm", "Could not locate MA0007.2.jaspar test file");
        } else {
            try {
                const PWM pwm = load_jaspar_pwm(jaspar_path);
                bool valid = (pwm.width > 0);
                for (size_t pos = 0; pos < pwm.width && valid; ++pos) {
                    double col_sum = 0.0;
                    for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
                        col_sum += pwm.at(b, pos);
                    }
                    if (!approx_equal(col_sum, 1.0, 1e-4)) {
                        valid = false;
                    }
                }

                if (valid) {
                    test_pass("load_jaspar_pwm loads normalized PWM");
                } else {
                    test_fail("load_jaspar_pwm normalization",
                              "Loaded PWM has invalid width or non-normalized columns");
                }
            } catch (const std::exception& e) {
                test_fail("load_jaspar_pwm exception", e.what());
            }
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
    test_markov_background();
    test_seed_discovery_utils();
    test_em_refinement_and_masking();
    test_io_utils();
    
    std::cout << "\n=============================================" << std::endl;
    std::cout << "  Results: " << tests_passed << " passed, " 
              << tests_failed << " failed" << std::endl;
    std::cout << "=============================================" << std::endl;
    
    return (tests_failed > 0) ? 1 : 0;
}
