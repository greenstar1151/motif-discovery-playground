/**
 * @file pwm_scoring.cpp
 * @brief Implementation of PWM Scoring (Core Operation 4)
 */

#include "motif/pwm_scoring.hpp"
#include <cmath>
#include <limits>

namespace motif {

double score_kmer(
    const std::string& kmer,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    // Check length match
    if (kmer.length() != pwm.width) {
        return -std::numeric_limits<double>::infinity();
    }
    
    double score = 0.0;
    
    for (size_t pos = 0; pos < pwm.width; ++pos) {
        size_t base_idx = base_to_index(kmer[pos]);
        
        // Invalid base (e.g., 'N')
        if (base_idx == SIZE_MAX) {
            return -std::numeric_limits<double>::infinity();
        }
        
        double pwm_prob = pwm.matrix[base_idx][pos];
        double bg_prob = background[base_idx];
        
        // Avoid log(0)
        if (pwm_prob <= 0.0 || bg_prob <= 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
        
        // Log-odds score
        score += std::log2(pwm_prob / bg_prob);
    }
    
    return score;
}

std::vector<double> scan_sequence(
    const Sequence& sequence,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    std::vector<double> scores;
    
    if (sequence.length() < pwm.width) {
        return scores;
    }
    
    size_t num_positions = sequence.length() - pwm.width + 1;
    scores.reserve(num_positions);
    
    for (size_t i = 0; i < num_positions; ++i) {
        std::string kmer = sequence.substr(i, pwm.width);
        scores.push_back(score_kmer(kmer, pwm, background));
    }
    
    return scores;
}

std::pair<int, double> find_best_match(
    const Sequence& sequence,
    const PWM& pwm,
    const std::array<double, ALPHABET_SIZE>& background
) {
    if (sequence.length() < pwm.width) {
        return {-1, -std::numeric_limits<double>::infinity()};
    }
    
    int best_pos = -1;
    double best_score = -std::numeric_limits<double>::infinity();
    
    size_t num_positions = sequence.length() - pwm.width + 1;
    
    for (size_t i = 0; i < num_positions; ++i) {
        std::string kmer = sequence.substr(i, pwm.width);
        double score = score_kmer(kmer, pwm, background);
        
        if (score > best_score) {
            best_score = score;
            best_pos = static_cast<int>(i);
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
    int hit_count = 0;
    
    for (const auto& seq : sequences) {
        auto [pos, score] = find_best_match(seq, pwm, background);
        if (score >= threshold) {
            hit_count++;
        }
    }
    
    return hit_count;
}

} // namespace motif
