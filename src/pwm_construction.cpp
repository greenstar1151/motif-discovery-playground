/**
 * @file pwm_construction.cpp
 * @brief Implementation of PWM Construction (Core Operation 3)
 */

#include "motif/pwm_construction.hpp"
#include <iostream>
#include <iomanip>

namespace motif {

PWM build_pwm_from_seed(
    const std::string& seed,
    const SequenceList& sequences,
    double pseudocount
) {
    // Collect all sites where seed is found
    std::vector<std::string> aligned_sites;
    
    for (const auto& seq : sequences) {
        size_t pos = seq.find(seed);
        if (pos != std::string::npos) {
            aligned_sites.push_back(seq.substr(pos, seed.length()));
        }
    }
    
    return build_pwm_from_sites(aligned_sites, pseudocount);
}

PWM build_pwm_from_sites(
    const std::vector<std::string>& aligned_sites,
    double pseudocount
) {
    PWM pwm;
    
    if (aligned_sites.empty()) {
        return pwm;
    }
    
    // Determine width from first site
    size_t width = aligned_sites[0].length();
    if (width > MAX_MOTIF_WIDTH) {
        width = MAX_MOTIF_WIDTH;
    }
    pwm.width = width;
    
    // Initialize frequency matrix with pseudocounts
    std::array<std::array<double, MAX_MOTIF_WIDTH>, ALPHABET_SIZE> freq = {};
    for (size_t base = 0; base < ALPHABET_SIZE; ++base) {
        for (size_t pos = 0; pos < width; ++pos) {
            freq[base][pos] = pseudocount;
        }
    }
    
    // Count base frequencies at each position
    size_t valid_sites = 0;
    for (const auto& site : aligned_sites) {
        if (site.length() < width) {
            continue;
        }
        
        bool valid = true;
        for (size_t pos = 0; pos < width; ++pos) {
            size_t base_idx = base_to_index(site[pos]);
            if (base_idx == SIZE_MAX) {
                valid = false;
                break;
            }
        }
        
        if (valid) {
            valid_sites++;
            for (size_t pos = 0; pos < width; ++pos) {
                size_t base_idx = base_to_index(site[pos]);
                freq[base_idx][pos] += 1.0;
            }
        }
    }
    
    // Normalize to probabilities
    double total = static_cast<double>(valid_sites) + pseudocount * ALPHABET_SIZE;
    
    for (size_t pos = 0; pos < width; ++pos) {
        for (size_t base = 0; base < ALPHABET_SIZE; ++base) {
            pwm.matrix[base][pos] = freq[base][pos] / total;
        }
    }
    
    return pwm;
}

std::string get_consensus(const PWM& pwm) {
    std::string consensus;
    consensus.reserve(pwm.width);
    
    for (size_t pos = 0; pos < pwm.width; ++pos) {
        size_t best_base = 0;
        double best_prob = pwm.matrix[0][pos];
        
        for (size_t base = 1; base < ALPHABET_SIZE; ++base) {
            if (pwm.matrix[base][pos] > best_prob) {
                best_prob = pwm.matrix[base][pos];
                best_base = base;
            }
        }
        
        consensus += index_to_base(best_base);
    }
    
    return consensus;
}

void print_pwm(const PWM& pwm, int precision) {
    if (pwm.width == 0) {
        std::cout << "  (Empty PWM)" << std::endl;
        return;
    }
    
    // Print header
    std::cout << "      ";
    for (size_t pos = 0; pos < pwm.width; ++pos) {
        std::cout << std::setw(precision + 3) << ("P" + std::to_string(pos + 1));
    }
    std::cout << std::endl;
    
    // Print each base row
    for (size_t base = 0; base < ALPHABET_SIZE; ++base) {
        std::cout << "  " << index_to_base(base) << ": ";
        for (size_t pos = 0; pos < pwm.width; ++pos) {
            std::cout << std::fixed << std::setprecision(precision) 
                      << std::setw(precision + 3) << pwm.matrix[base][pos];
        }
        std::cout << std::endl;
    }
}

} // namespace motif
