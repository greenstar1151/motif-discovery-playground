/**
 * @file pwm_construction.cpp
 * @brief Implementation of PWM construction.
 */

#include "motif/pwm_construction.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>

namespace motif {

namespace {

size_t hamming_distance(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) {
        return std::numeric_limits<size_t>::max();
    }
    size_t dist = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            ++dist;
        }
    }
    return dist;
}

} // namespace

PWM build_pwm_from_seed(
    const std::string& seed,
    const SequenceList& sequences,
    double pseudocount
) {
    std::vector<std::string> aligned_sites;
    if (seed.empty()) {
        return {};
    }

    for (const auto& seq : sequences) {
        if (seq.size() < seed.size()) {
            continue;
        }

        size_t best_hd = std::numeric_limits<size_t>::max();
        size_t best_pos = std::numeric_limits<size_t>::max();

        for (size_t pos = 0; pos + seed.size() <= seq.size(); ++pos) {
            const std::string window = seq.substr(pos, seed.size());
            if (!is_valid_kmer(window)) {
                continue;
            }

            const size_t hd = hamming_distance(window, seed);
            if (hd < best_hd) {
                best_hd = hd;
                best_pos = pos;
                if (hd == 0) {
                    break;
                }
            }
        }

        if (best_pos != std::numeric_limits<size_t>::max() && best_hd <= 1) {
            aligned_sites.push_back(seq.substr(best_pos, seed.size()));
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

    const double safe_pseudocount = (pseudocount < 0.0) ? 0.0 : pseudocount;
    size_t width = aligned_sites.front().size();
    width = std::min(width, MAX_MOTIF_WIDTH);
    if (width == 0) {
        return pwm;
    }
    pwm.width = width;

    std::array<std::array<double, MAX_MOTIF_WIDTH>, ALPHABET_SIZE> counts{};
    std::array<double, MAX_MOTIF_WIDTH> column_totals{};

    for (size_t pos = 0; pos < width; ++pos) {
        for (size_t base = 0; base < ALPHABET_SIZE; ++base) {
            counts[base][pos] = safe_pseudocount;
        }
        column_totals[pos] = safe_pseudocount * static_cast<double>(ALPHABET_SIZE);
    }

    for (const auto& site : aligned_sites) {
        if (site.size() < width) {
            continue;
        }

        bool valid = true;
        for (size_t pos = 0; pos < width; ++pos) {
            if (base_to_index(site[pos]) == SIZE_MAX) {
                valid = false;
                break;
            }
        }
        if (!valid) {
            continue;
        }

        for (size_t pos = 0; pos < width; ++pos) {
            const size_t base_idx = base_to_index(site[pos]);
            counts[base_idx][pos] += 1.0;
            column_totals[pos] += 1.0;
        }
    }

    for (size_t pos = 0; pos < width; ++pos) {
        const double denom = (column_totals[pos] > 0.0) ? column_totals[pos] : 1.0;
        for (size_t base = 0; base < ALPHABET_SIZE; ++base) {
            pwm.matrix[base][pos] = counts[base][pos] / denom;
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
        consensus.push_back(index_to_base(best_base));
    }

    return consensus;
}

void print_pwm(const PWM& pwm, int precision) {
    if (pwm.width == 0) {
        std::cout << "  (Empty PWM)\n";
        return;
    }

    std::cout << "      ";
    for (size_t pos = 0; pos < pwm.width; ++pos) {
        std::cout << std::setw(precision + 4) << ("P" + std::to_string(pos + 1));
    }
    std::cout << '\n';

    for (size_t base = 0; base < ALPHABET_SIZE; ++base) {
        std::cout << "  " << index_to_base(base) << ": ";
        for (size_t pos = 0; pos < pwm.width; ++pos) {
            std::cout << std::fixed << std::setprecision(precision)
                      << std::setw(precision + 4) << pwm.matrix[base][pos];
        }
        std::cout << '\n';
    }
}

} // namespace motif
