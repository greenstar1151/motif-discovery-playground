/**
 * @file markov.hpp
 * @brief 3rd-order Markov background model for DNA sequences.
 *
 * Trains on control sequences and provides conditional probabilities
 * P(base | 3-base context) used as the background model for LLR scoring.
 */

#ifndef MOTIF_MARKOV_HPP
#define MOTIF_MARKOV_HPP

#include "motif/types.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace motif {

/// Minimum probability to avoid log(0).
constexpr double kMinProb = 1e-12;

/**
 * @class MarkovOrder3
 * @brief Order-3 Markov background model for DNA.
 *
 * Context is the preceding 3 bases; falls back to marginal
 * base frequencies for positions 0–2.
 */
class MarkovOrder3 {
public:
    /**
     * @brief Train the model from a set of sequences.
     * @param sequences Training sequences (typically control set).
     * @param pseudocount Laplace pseudocount per cell (default 1.0).
     */
    void train(const SequenceList& sequences, double pseudocount = 1.0);

    /**
     * @brief Conditional probability P(seq[pos] | seq[pos-3..pos-1]).
     * @return Probability in (0, 1]. Returns 0.25 for invalid bases.
     */
    double conditional_prob(const Sequence& seq, size_t pos) const;

private:
    std::array<std::array<double, ALPHABET_SIZE>, 64> transition_counts_{};
    std::array<std::array<double, ALPHABET_SIZE>, 64> transition_probs_{};
    std::array<double, ALPHABET_SIZE> base_probs_ = {0.25, 0.25, 0.25, 0.25};
};

} // namespace motif

#endif // MOTIF_MARKOV_HPP
