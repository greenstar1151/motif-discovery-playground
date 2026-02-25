/**
 * @file pwm_scoring.cpp
 * @brief Implementation of PWM scoring utilities.
 */

#include "motif/pwm_scoring.hpp"
#include <cmath>
#include <limits>

namespace motif {

namespace {

constexpr double kMinProb = 1e-12;

double score_window(
    const Sequence& sequence,
    size_t start,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    double score = 0.0;
    for (size_t pos = 0; pos < pwm.width; ++pos) {
        const size_t base_idx = base_to_index(sequence[start + pos]);
        if (base_idx == SIZE_MAX) {
            return -std::numeric_limits<double>::infinity();
        }

        const double pwm_prob = std::max(pwm.matrix[base_idx][pos], kMinProb);
        const double bg_prob = std::max(background[base_idx], kMinProb);
        score += std::log2(pwm_prob / bg_prob);
    }
    return score;
}

} // namespace

double score_kmer(
    const std::string& kmer,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    if (kmer.size() != pwm.width || pwm.width == 0) {
        return -std::numeric_limits<double>::infinity();
    }
    return score_window(kmer, 0, pwm, background);
}

std::vector<double> scan_sequence(
    const Sequence& sequence,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    std::vector<double> scores;
    if (pwm.width == 0 || sequence.size() < pwm.width) {
        return scores;
    }

    const size_t n_positions = sequence.size() - pwm.width + 1;
    scores.reserve(n_positions);
    for (size_t pos = 0; pos < n_positions; ++pos) {
        scores.push_back(score_window(sequence, pos, pwm, background));
    }
    return scores;
}

std::pair<int, double> find_best_match(
    const Sequence& sequence,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    if (pwm.width == 0 || sequence.size() < pwm.width) {
        return {-1, -std::numeric_limits<double>::infinity()};
    }

    int best_pos = -1;
    double best_score = -std::numeric_limits<double>::infinity();
    const size_t n_positions = sequence.size() - pwm.width + 1;

    for (size_t pos = 0; pos < n_positions; ++pos) {
        const double score = score_window(sequence, pos, pwm, background);
        if (score > best_score) {
            best_score = score;
            best_pos = static_cast<int>(pos);
        }
    }

    return {best_pos, best_score};
}

int count_hits(
    const SequenceList& sequences,
    const PWM& pwm,
    double threshold,
    const std::array<double, ALPHABET_SIZE>& background
) {
    int hits = 0;
    for (const auto& seq : sequences) {
        const auto [_, best_score] = find_best_match(seq, pwm, background);
        if (best_score >= threshold) {
            ++hits;
        }
    }
    return hits;
}

} // namespace motif
