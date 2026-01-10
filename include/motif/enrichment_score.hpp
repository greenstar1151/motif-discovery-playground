/**
 * @file enrichment_score.hpp
 * @brief Enrichment Score Calculation
 * 
 * Core Operation 2: Compare k-mer frequencies between primary and control
 * sets to identify statistically significant patterns.
 * 
 * Uses simple ratio difference for clarity. More sophisticated tests
 * (Fisher's exact, binomial) can be added for production use.
 * 
 * Complexity: O(4^k) in worst case, typically much smaller
 */

#ifndef MOTIF_ENRICHMENT_SCORE_HPP
#define MOTIF_ENRICHMENT_SCORE_HPP

#include "motif/types.hpp"

namespace motif {

/**
 * @brief Calculate enrichment scores for all k-mers
 * 
 * Score = (primary_rate) - (control_rate)
 * where rate = count / total_sequences
 * 
 * Positive score indicates enrichment in primary set.
 * 
 * @param pos_counts K-mer counts in primary (positive) set
 * @param neg_counts K-mer counts in control (negative) set
 * @param pos_total Total sequences in primary set
 * @param neg_total Total sequences in control set
 * @return Map of k-mer -> enrichment score
 * 
 * @example
 *   If "TATATA" appears in 90/200 primary and 5/200 control:
 *   score = (90/200) - (5/200) = 0.45 - 0.025 = 0.425
 */
EnrichmentScores calculate_enrichment_scores(
    const KmerCounts& pos_counts,
    const KmerCounts& neg_counts,
    int pos_total,
    int neg_total
);

/**
 * @brief Find the k-mer with highest enrichment score
 * 
 * @param scores Enrichment scores for all k-mers
 * @return Pair of (best_kmer, best_score), or ("", 0.0) if empty
 */
std::pair<std::string, double> find_best_seed(const EnrichmentScores& scores);

/**
 * @brief Get top N k-mers by enrichment score
 * 
 * @param scores Enrichment scores for all k-mers
 * @param n Number of top k-mers to return
 * @return Vector of (kmer, score) pairs sorted by score descending
 */
std::vector<std::pair<std::string, double>> get_top_seeds(
    const EnrichmentScores& scores, 
    size_t n
);

} // namespace motif

#endif // MOTIF_ENRICHMENT_SCORE_HPP
