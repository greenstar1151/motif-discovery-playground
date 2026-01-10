/**
 * @file pwm_construction.hpp
 * @brief Position Weight Matrix Construction
 * 
 * Core Operation 3: Build a PWM from aligned sequence sites.
 * The PWM captures the positional base preferences of a motif.
 * 
 * Complexity: O(M × W) where M = number of sites, W = motif width
 */

#ifndef MOTIF_PWM_CONSTRUCTION_HPP
#define MOTIF_PWM_CONSTRUCTION_HPP

#include "motif/types.hpp"

namespace motif {

/**
 * @brief Build PWM from a seed k-mer and sequences
 * 
 * Finds all exact matches of the seed in sequences and builds
 * a Position Frequency Matrix, then normalizes to probabilities.
 * 
 * @param seed The consensus k-mer to use as anchor
 * @param sequences DNA sequences to search
 * @param pseudocount Laplace smoothing constant (default 0.1)
 * @return PWM with normalized probabilities
 * 
 * @example
 *   seed = "TATATA"
 *   If found in 100 sequences, the PWM will show high probability
 *   for T at positions 1,3,5 and A at positions 2,4,6.
 */
PWM build_pwm_from_seed(
    const std::string& seed,
    const SequenceList& sequences,
    double pseudocount = 0.1
);

/**
 * @brief Build PWM from aligned sites directly
 * 
 * @param aligned_sites Vector of equal-length DNA strings
 * @param pseudocount Laplace smoothing constant (default 0.1)
 * @return PWM with normalized probabilities
 */
PWM build_pwm_from_sites(
    const std::vector<std::string>& aligned_sites,
    double pseudocount = 0.1
);

/**
 * @brief Get consensus sequence from PWM
 * 
 * For each position, selects the base with highest probability.
 * 
 * @param pwm Position Weight Matrix
 * @return Consensus sequence string
 */
std::string get_consensus(const PWM& pwm);

/**
 * @brief Print PWM in human-readable format
 * 
 * @param pwm Position Weight Matrix
 * @param precision Decimal places for probabilities (default 3)
 */
void print_pwm(const PWM& pwm, int precision = 3);

} // namespace motif

#endif // MOTIF_PWM_CONSTRUCTION_HPP
