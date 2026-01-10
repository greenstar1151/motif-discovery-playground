/**
 * @file sequence_masking.hpp
 * @brief Sequence Masking Operation
 * 
 * Core Operation 5: Remove discovered motifs from sequences by replacing
 * with 'N' characters. This allows iterative discovery of multiple motifs.
 * 
 * Complexity: O(N × L) where N = number of sequences, L = average length
 */

#ifndef MOTIF_SEQUENCE_MASKING_HPP
#define MOTIF_SEQUENCE_MASKING_HPP

#include "motif/types.hpp"

namespace motif {

/**
 * @brief Mask all occurrences of a pattern in sequences
 * 
 * Replaces exact matches of the pattern with 'N' characters.
 * 
 * @param sequences List of DNA sequences
 * @param pattern Pattern to mask (exact string match)
 * @return New sequence list with pattern occurrences masked
 * 
 * @example
 *   sequence = "ACGTATATACGT"
 *   pattern = "TATATA"
 *   result = "ACGNNNNNNNCGT"
 */
SequenceList mask_pattern(const SequenceList& sequences, const std::string& pattern);

/**
 * @brief Mask in place (modifies input sequences)
 * 
 * @param sequences List of DNA sequences (modified in place)
 * @param pattern Pattern to mask
 * @return Number of occurrences masked
 */
int mask_pattern_inplace(SequenceList& sequences, const std::string& pattern);

/**
 * @brief Mask positions based on PWM score threshold
 * 
 * Finds positions scoring above threshold and masks them.
 * Useful for masking fuzzy matches.
 * 
 * @param sequences List of DNA sequences
 * @param pwm Position Weight Matrix
 * @param threshold Minimum score to mask
 * @return New sequence list with high-scoring positions masked
 */
SequenceList mask_by_pwm(
    const SequenceList& sequences,
    const PWM& pwm,
    double threshold
);

/**
 * @brief Count masked positions in a sequence
 * 
 * @param sequence DNA sequence
 * @return Number of 'N' characters
 */
size_t count_masked(const Sequence& sequence);

/**
 * @brief Calculate fraction of sequence that is masked
 * 
 * @param sequences List of DNA sequences
 * @return Fraction of total bases that are 'N'
 */
double masked_fraction(const SequenceList& sequences);

} // namespace motif

#endif // MOTIF_SEQUENCE_MASKING_HPP
