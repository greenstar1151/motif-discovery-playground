/**
 * @file pwm_scoring.hpp
 * @brief PWM Scoring Operation
 * 
 * Core Operation 4: Score k-mers or sequence positions against a PWM.
 * Uses log-odds scoring against a background model.
 * 
 * Complexity: O(W) per k-mer where W = motif width
 */

#ifndef MOTIF_PWM_SCORING_HPP
#define MOTIF_PWM_SCORING_HPP

#include "motif/types.hpp"

namespace motif {

/**
 * @brief Default uniform background model (0.25 for each base)
 */
constexpr std::array<double, ALPHABET_SIZE> UNIFORM_BACKGROUND = {0.25, 0.25, 0.25, 0.25};

/**
 * @brief Score a k-mer against a PWM using log-odds
 * 
 * Score = Σ log2(PWM[base][pos] / background[base])
 * 
 * Higher scores indicate better match to the motif model.
 * 
 * @param kmer K-mer string to score (must match PWM width)
 * @param pwm Position Weight Matrix
 * @param background Background base probabilities (default uniform)
 * @return Log-odds score
 * 
 * @example
 *   If PWM strongly prefers "TATATA":
 *   score("TATATA") >> score("ACGTAC")
 */
double score_kmer(
    const std::string& kmer,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background = UNIFORM_BACKGROUND
);

/**
 * @brief Scan a sequence and return scores at all positions
 * 
 * @param sequence DNA sequence to scan
 * @param pwm Position Weight Matrix
 * @param background Background base probabilities
 * @return Vector of scores at each valid position
 */
std::vector<double> scan_sequence(
    const Sequence& sequence,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background = UNIFORM_BACKGROUND
);

/**
 * @brief Find the best scoring position in a sequence
 * 
 * @param sequence DNA sequence to scan
 * @param pwm Position Weight Matrix
 * @param background Background base probabilities
 * @return Pair of (best_position, best_score), or (-1, -inf) if no valid positions
 */
std::pair<int, double> find_best_match(
    const Sequence& sequence,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background = UNIFORM_BACKGROUND
);

/**
 * @brief Count sequences with a score above threshold
 * 
 * @param sequences List of DNA sequences
 * @param pwm Position Weight Matrix
 * @param threshold Minimum score to count
 * @param background Background base probabilities
 * @return Number of sequences with at least one position above threshold
 */
int count_hits(
    const SequenceList& sequences,
    const PWM& pwm,
    double threshold,
    const std::array<double, ALPHABET_SIZE>& background = UNIFORM_BACKGROUND
);

} // namespace motif

#endif // MOTIF_PWM_SCORING_HPP
