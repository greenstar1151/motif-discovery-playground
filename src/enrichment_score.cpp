/**
 * @file enrichment_score.cpp
 * @brief Implementation of enrichment-scoring helpers.
 */

#include "motif/enrichment_score.hpp"
#include <algorithm>
#include <set>

namespace motif {

EnrichmentScores calculate_enrichment_scores(
    const KmerCounts& pos_counts,
    const KmerCounts& neg_counts,
    int pos_total,
    int neg_total
) {
    EnrichmentScores scores;

    if (pos_total <= 0 || neg_total <= 0) {
        return scores;
    }

    std::set<std::string> all_kmers;
    for (const auto& [kmer, _] : pos_counts) {
        all_kmers.insert(kmer);
    }
    for (const auto& [kmer, _] : neg_counts) {
        all_kmers.insert(kmer);
    }

    const double inv_pos_total = 1.0 / static_cast<double>(pos_total);
    const double inv_neg_total = 1.0 / static_cast<double>(neg_total);

    for (const auto& kmer : all_kmers) {
        const auto pos_it = pos_counts.find(kmer);
        const auto neg_it = neg_counts.find(kmer);

        const int pos_count = (pos_it == pos_counts.end()) ? 0 : pos_it->second;
        const int neg_count = (neg_it == neg_counts.end()) ? 0 : neg_it->second;

        const double pos_rate = static_cast<double>(pos_count) * inv_pos_total;
        const double neg_rate = static_cast<double>(neg_count) * inv_neg_total;
        scores[kmer] = pos_rate - neg_rate;
    }

    return scores;
}

std::pair<std::string, double> find_best_seed(const EnrichmentScores& scores) {
    if (scores.empty()) {
        return {"", 0.0};
    }

    auto best_it = std::max_element(
        scores.begin(),
        scores.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; }
    );
    return {best_it->first, best_it->second};
}

std::vector<std::pair<std::string, double>> get_top_seeds(
    const EnrichmentScores& scores,
    size_t n
) {
    std::vector<std::pair<std::string, double>> sorted(scores.begin(), scores.end());
    std::sort(
        sorted.begin(),
        sorted.end(),
        [](const auto& left, const auto& right) {
            if (left.second != right.second) {
                return left.second > right.second;
            }
            return left.first < right.first;
        }
    );

    if (sorted.size() > n) {
        sorted.resize(n);
    }
    return sorted;
}

} // namespace motif
