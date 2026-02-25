/**
 * @file sequence_masking.cpp
 * @brief Implementation of sequence-masking helpers.
 */

#include "motif/sequence_masking.hpp"
#include "motif/pwm_scoring.hpp"

namespace motif {

SequenceList mask_pattern(const SequenceList& sequences, const std::string& pattern) {
    if (pattern.empty()) {
        return sequences;
    }

    SequenceList result;
    result.reserve(sequences.size());
    const std::string mask(pattern.size(), MASK_CHAR);

    for (const auto& seq : sequences) {
        std::string masked = seq;
        size_t pos = 0;
        while ((pos = masked.find(pattern, pos)) != std::string::npos) {
            masked.replace(pos, pattern.size(), mask);
            pos += pattern.size();
        }
        result.push_back(std::move(masked));
    }

    return result;
}

int mask_pattern_inplace(SequenceList& sequences, const std::string& pattern) {
    if (pattern.empty()) {
        return 0;
    }

    int replaced = 0;
    const std::string mask(pattern.size(), MASK_CHAR);

    for (auto& seq : sequences) {
        size_t pos = 0;
        while ((pos = seq.find(pattern, pos)) != std::string::npos) {
            seq.replace(pos, pattern.size(), mask);
            pos += pattern.size();
            ++replaced;
        }
    }

    return replaced;
}

SequenceList mask_by_pwm(
    const SequenceList& sequences,
    const PWM& pwm,
    double threshold
) {
    if (pwm.width == 0) {
        return sequences;
    }

    SequenceList result;
    result.reserve(sequences.size());
    const std::string mask(pwm.width, MASK_CHAR);

    for (const auto& seq : sequences) {
        std::string masked = seq;

        if (seq.size() >= pwm.width) {
            std::vector<size_t> to_mask;
            to_mask.reserve(seq.size() - pwm.width + 1);

            for (size_t pos = 0; pos + pwm.width <= seq.size(); ++pos) {
                const std::string window = seq.substr(pos, pwm.width);
                const double score = score_kmer(window, pwm);
                if (score >= threshold) {
                    to_mask.push_back(pos);
                }
            }

            for (auto it = to_mask.rbegin(); it != to_mask.rend(); ++it) {
                masked.replace(*it, pwm.width, mask);
            }
        }

        result.push_back(std::move(masked));
    }

    return result;
}

size_t count_masked(const Sequence& sequence) {
    size_t count = 0;
    for (char base : sequence) {
        if (base == MASK_CHAR) {
            ++count;
        }
    }
    return count;
}

double masked_fraction(const SequenceList& sequences) {
    size_t total = 0;
    size_t masked = 0;
    for (const auto& seq : sequences) {
        total += seq.size();
        masked += count_masked(seq);
    }
    return (total == 0) ? 0.0 : static_cast<double>(masked) / static_cast<double>(total);
}

} // namespace motif
