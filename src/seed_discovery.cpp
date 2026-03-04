/**
 * @file seed_discovery.cpp
 * @brief Implementation of multi-k z-test seed discovery.
 */

#include "motif/seed_discovery.hpp"
#include "motif/kmer_counting.hpp"

#include <algorithm>
#include <cmath>
#include <set>

namespace motif {

// ─────────────────────────────────────────────────────────────────────────────
// z_score_two_proportion
// ─────────────────────────────────────────────────────────────────────────────

double z_score_two_proportion(int x1, int n1, int x2, int n2) {
    if (n1 <= 0 || n2 <= 0) {
        return 0.0;
    }

    const double p1 = static_cast<double>(x1) / static_cast<double>(n1);
    const double p2 = static_cast<double>(x2) / static_cast<double>(n2);
    const double pooled =
        static_cast<double>(x1 + x2) / static_cast<double>(n1 + n2);
    const double var =
        pooled * (1.0 - pooled) * (1.0 / n1 + 1.0 / n2);

    if (var <= 0.0) {
        return 0.0;
    }
    return (p1 - p2) / std::sqrt(var);
}

// ─────────────────────────────────────────────────────────────────────────────
// hamming_distance
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// collect_hd_sites
// ─────────────────────────────────────────────────────────────────────────────

std::vector<std::string> collect_hd_sites(
    const std::string& seed,
    const SequenceList& sequences,
    size_t target_width,
    size_t max_hd
) {
    std::vector<std::string> sites;
    if (seed.empty()) {
        return sites;
    }

    const size_t width = (target_width > 0) ? target_width : seed.size();
    if (width == 0 || width > MAX_MOTIF_WIDTH) {
        return sites;
    }

    sites.reserve(sequences.size());

    for (const auto& seq : sequences) {
        if (seq.size() < seed.size() || seq.size() < width) {
            continue;
        }

        size_t best_hd = std::numeric_limits<size_t>::max();
        size_t best_seed_pos = std::numeric_limits<size_t>::max();

        for (size_t pos = 0; pos + seed.size() <= seq.size(); ++pos) {
            const std::string window = seq.substr(pos, seed.size());
            if (!is_valid_kmer(window)) {
                continue;
            }

            const size_t hd = hamming_distance(window, seed);
            if (hd < best_hd) {
                best_hd = hd;
                best_seed_pos = pos;
                if (hd == 0) {
                    break;
                }
            }
        }

        if (best_seed_pos == std::numeric_limits<size_t>::max() || best_hd > max_hd) {
            continue;
        }

        // Center width around the seed match when width != seed.size()
        size_t start = best_seed_pos;
        if (width != seed.size()) {
            const size_t center = best_seed_pos + seed.size() / 2;
            if (center >= width / 2) {
                start = center - width / 2;
            } else {
                start = 0;
            }
            if (start + width > seq.size()) {
                start = seq.size() - width;
            }
        }

        const std::string site = seq.substr(start, width);
        if (is_valid_kmer(site)) {
            sites.push_back(site);
        }
    }

    return sites;
}

// ─────────────────────────────────────────────────────────────────────────────
// build_seed_pool
// ─────────────────────────────────────────────────────────────────────────────

std::vector<SeedCandidate> build_seed_pool(
    const SequenceList& primary,
    const SequenceList& control,
    const std::vector<size_t>& k_values,
    double z_sub_threshold,
    size_t top_per_k
) {
    std::vector<SeedCandidate> pool;

    for (size_t k : k_values) {
        const KmerCounts pos_counts = count_kmers(primary, k);
        const KmerCounts neg_counts = count_kmers(control, k);

        if (pos_counts.empty() && neg_counts.empty()) {
            continue;
        }

        std::set<std::string> all_kmers;
        for (const auto& [kmer, _] : pos_counts) {
            all_kmers.insert(kmer);
        }
        for (const auto& [kmer, _] : neg_counts) {
            all_kmers.insert(kmer);
        }

        std::vector<SeedCandidate> local;
        local.reserve(all_kmers.size());

        for (const auto& kmer : all_kmers) {
            const int x1 = pos_counts.count(kmer) ? pos_counts.at(kmer) : 0;
            const int x2 = neg_counts.count(kmer) ? neg_counts.at(kmer) : 0;
            const double z = z_score_two_proportion(
                x1,
                static_cast<int>(primary.size()),
                x2,
                static_cast<int>(control.size())
            );

            if (z >= z_sub_threshold) {
                local.push_back({kmer, k, z, x1, x2});
            }
        }

        std::sort(
            local.begin(), local.end(),
            [](const SeedCandidate& a, const SeedCandidate& b) {
                if (a.z_score != b.z_score) {
                    return a.z_score > b.z_score;
                }
                return a.kmer < b.kmer;
            }
        );

        if (local.size() > top_per_k) {
            local.resize(top_per_k);
        }

        pool.insert(pool.end(), local.begin(), local.end());
    }

    // Global sort: z descending
    std::sort(
        pool.begin(), pool.end(),
        [](const SeedCandidate& a, const SeedCandidate& b) {
            if (a.z_score != b.z_score) {
                return a.z_score > b.z_score;
            }
            return a.kmer < b.kmer;
        }
    );

    return pool;
}

} // namespace motif
