/**
 * @file markov.cpp
 * @brief Implementation of MarkovOrder3 background model.
 */

#include "motif/markov.hpp"

namespace motif {

void MarkovOrder3::train(const SequenceList& sequences, double pseudocount) {
    std::array<double, ALPHABET_SIZE> base_counts{};
    base_counts.fill(pseudocount);

    for (auto& row : transition_counts_) {
        row.fill(pseudocount);
    }

    for (const auto& seq : sequences) {
        for (size_t pos = 0; pos < seq.size(); ++pos) {
            const size_t base_idx = base_to_index(seq[pos]);
            if (base_idx == SIZE_MAX) {
                continue;
            }

            base_counts[base_idx] += 1.0;

            if (pos < 3) {
                continue;
            }

            const size_t b1 = base_to_index(seq[pos - 3]);
            const size_t b2 = base_to_index(seq[pos - 2]);
            const size_t b3 = base_to_index(seq[pos - 1]);
            if (b1 == SIZE_MAX || b2 == SIZE_MAX || b3 == SIZE_MAX) {
                continue;
            }

            const size_t context = ((b1 * 4) + b2) * 4 + b3;
            transition_counts_[context][base_idx] += 1.0;
        }
    }

    double base_total = 0.0;
    for (double v : base_counts) {
        base_total += v;
    }
    for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
        base_probs_[b] = base_counts[b] / std::max(base_total, kMinProb);
    }

    for (size_t ctx = 0; ctx < 64; ++ctx) {
        double row_total = 0.0;
        for (double v : transition_counts_[ctx]) {
            row_total += v;
        }
        for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
            transition_probs_[ctx][b] =
                transition_counts_[ctx][b] / std::max(row_total, kMinProb);
        }
    }
}

double MarkovOrder3::conditional_prob(const Sequence& seq, size_t pos) const {
    if (pos >= seq.size()) {
        return 0.25;
    }

    const size_t base_idx = base_to_index(seq[pos]);
    if (base_idx == SIZE_MAX) {
        return 0.25;
    }

    if (pos < 3) {
        return std::max(base_probs_[base_idx], kMinProb);
    }

    const size_t b1 = base_to_index(seq[pos - 3]);
    const size_t b2 = base_to_index(seq[pos - 2]);
    const size_t b3 = base_to_index(seq[pos - 1]);
    if (b1 == SIZE_MAX || b2 == SIZE_MAX || b3 == SIZE_MAX) {
        return std::max(base_probs_[base_idx], kMinProb);
    }

    const size_t context = ((b1 * 4) + b2) * 4 + b3;
    return std::max(transition_probs_[context][base_idx], kMinProb);
}

} // namespace motif
