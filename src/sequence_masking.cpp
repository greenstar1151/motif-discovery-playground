/**
 * @file sequence_masking.cpp
 * @brief Implementation of Sequence Masking (Core Operation 5)
 */

#include "motif/sequence_masking.hpp"
#include "motif/pwm_scoring.hpp"

namespace motif {

SequenceList mask_pattern(const SequenceList& sequences, const std::string& pattern) {
    SequenceList result;
    result.reserve(sequences.size());
    
    std::string mask(pattern.length(), MASK_CHAR);
    
    for (const auto& seq : sequences) {
        std::string masked_seq = seq;
        
        // Replace all occurrences
        size_t pos = 0;
        while ((pos = masked_seq.find(pattern, pos)) != std::string::npos) {
            masked_seq.replace(pos, pattern.length(), mask);
            pos += pattern.length();  // Move past the replacement
        }
        
        result.push_back(masked_seq);
    }
    
    return result;
}

int mask_pattern_inplace(SequenceList& sequences, const std::string& pattern) {
    int total_masked = 0;
    std::string mask(pattern.length(), MASK_CHAR);
    
    for (auto& seq : sequences) {
        size_t pos = 0;
        while ((pos = seq.find(pattern, pos)) != std::string::npos) {
            seq.replace(pos, pattern.length(), mask);
            pos += pattern.length();
            total_masked++;
        }
    }
    
    return total_masked;
}

SequenceList mask_by_pwm(
    const SequenceList& sequences,
    const PWM& pwm,
    double threshold
) {
    SequenceList result;
    result.reserve(sequences.size());
    
    std::string mask(pwm.width, MASK_CHAR);
    
    for (const auto& seq : sequences) {
        std::string masked_seq = seq;
        
        if (seq.length() >= pwm.width) {
            // Scan sequence and mask positions above threshold
            // Process from end to start to avoid position shifting issues
            std::vector<size_t> positions_to_mask;
            
            for (size_t i = 0; i <= seq.length() - pwm.width; ++i) {
                std::string kmer = seq.substr(i, pwm.width);
                double score = score_kmer(kmer, pwm);
                
                if (score >= threshold) {
                    positions_to_mask.push_back(i);
                }
            }
            
            // Mask from end to start
            for (auto it = positions_to_mask.rbegin(); it != positions_to_mask.rend(); ++it) {
                masked_seq.replace(*it, pwm.width, mask);
            }
        }
        
        result.push_back(masked_seq);
    }
    
    return result;
}

size_t count_masked(const Sequence& sequence) {
    size_t count = 0;
    for (char c : sequence) {
        if (c == MASK_CHAR) {
            count++;
        }
    }
    return count;
}

double masked_fraction(const SequenceList& sequences) {
    size_t total_bases = 0;
    size_t masked_bases = 0;
    
    for (const auto& seq : sequences) {
        total_bases += seq.length();
        masked_bases += count_masked(seq);
    }
    
    if (total_bases == 0) {
        return 0.0;
    }
    
    return static_cast<double>(masked_bases) / total_bases;
}

} // namespace motif
