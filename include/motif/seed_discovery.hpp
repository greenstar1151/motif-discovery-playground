/**
 * @file seed_discovery.hpp
 * @brief Seed discovery via multi-k two-proportion z-test.
 *
 * Provides:
 *   - Two-proportion z-test scoring
 *   - Hamming distance computation
 *   - HD≤d site collection for PWM initialization
 *   - Multi-k seed pool construction
 */

#ifndef MOTIF_SEED_DISCOVERY_HPP
#define MOTIF_SEED_DISCOVERY_HPP

#include "motif/types.hpp"
#include <limits>
#include <string>
#include <vector>

namespace motif {

// ─────────────────────────────────────────────────────────────────────────────
// Data structures
// ─────────────────────────────────────────────────────────────────────────────

/** @brief A single seed candidate from z-test enrichment screening. */
struct SeedCandidate {
    std::string kmer;
    size_t width = 0;
    double z_score = 0.0;
    int pos_count = 0;
    int neg_count = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Free functions
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Two-proportion z-test score.
 *
 * z = (p1 - p2) / sqrt(pooled * (1 - pooled) * (1/n1 + 1/n2))
 */
double z_score_two_proportion(int x1, int n1, int x2, int n2);

/**
 * @brief Hamming distance between two equal-length strings.
 * @return Distance, or SIZE_MAX if lengths differ.
 */
size_t hamming_distance(const std::string& a, const std::string& b);

/**
 * @brief Collect best HD≤max_hd sites from sequences, one per sequence.
 *
 * For each sequence the best-matching (lowest HD) window to @p seed is
 * selected.  If @p target_width differs from seed length the window is
 * centered on the seed match and extended/trimmed to @p target_width.
 *
 * @param seed            Anchor k-mer.
 * @param sequences       Sequences to scan.
 * @param target_width    Desired site width (0 = same as seed).
 * @param max_hd          Maximum Hamming distance to accept.
 * @return Collected site strings.
 */
std::vector<std::string> collect_hd_sites(
    const std::string& seed,
    const SequenceList& sequences,
    size_t target_width = 0,
    size_t max_hd = 1
);

/**
 * @brief Build a ranked pool of seed candidates across multiple k values.
 *
 * For every k in @p k_values, k-mers are counted in both sets, scored
 * with a two-proportion z-test, filtered by @p z_sub_threshold, and the
 * top @p top_per_k are kept.  The final pool is sorted by z descending.
 *
 * @param primary          Primary (positive) sequence set.
 * @param control          Control (negative) sequence set.
 * @param k_values         Seed widths to scan.
 * @param z_sub_threshold  Minimum z-score to retain a candidate.
 * @param top_per_k        Maximum candidates to retain per k value.
 * @return Merged, sorted seed pool.
 */
std::vector<SeedCandidate> build_seed_pool(
    const SequenceList& primary,
    const SequenceList& control,
    const std::vector<size_t>& k_values,
    double z_sub_threshold,
    size_t top_per_k = 25
);

} // namespace motif

#endif // MOTIF_SEED_DISCOVERY_HPP
