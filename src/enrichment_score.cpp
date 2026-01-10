/**
 * @file enrichment_score.cpp
 * @brief Implementation of Enrichment Score Calculation (Core Operation 2)
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
    
    // Collect all unique k-mers from both sets
    std::set<std::string> all_kmers;
    for (const auto& [kmer, _] : pos_counts) {
        all_kmers.insert(kmer);
    }
    for (const auto& [kmer, _] : neg_counts) {
        all_kmers.insert(kmer);
    }
    
    // Calculate enrichment score for each k-mer
    for (const auto& kmer : all_kmers) {
        // Get counts (default to 0 if not present)
        int pos_count = 0;
        int neg_count = 0;
        
        auto pos_it = pos_counts.find(kmer);
        if (pos_it != pos_counts.end()) {
            pos_count = pos_it->second;
        }
        
        auto neg_it = neg_counts.find(kmer);
        if (neg_it != neg_counts.end()) {
            neg_count = neg_it->second;
        }
        
        // Calculate rates
        double pos_rate = static_cast<double>(pos_count) / pos_total;
        double neg_rate = static_cast<double>(neg_count) / neg_total;
        
        // Enrichment score = positive rate - negative rate
        scores[kmer] = pos_rate - neg_rate;
    }
    
    return scores;
}

std::pair<std::string, double> find_best_seed(const EnrichmentScores& scores) {
    if (scores.empty()) {
        return {"", 0.0};
    }
    
    std::string best_kmer;
    double best_score = -1e9;  // Very negative initial value
    
    for (const auto& [kmer, score] : scores) {
        if (score > best_score) {
            best_score = score;
            best_kmer = kmer;
        }
    }
    
    return {best_kmer, best_score};
}

std::vector<std::pair<std::string, double>> get_top_seeds(
    const EnrichmentScores& scores, 
    size_t n
) {
    // Convert to vector for sorting
    std::vector<std::pair<std::string, double>> sorted_scores(
        scores.begin(), 
        scores.end()
    );
    
    // Sort by score descending
    std::sort(sorted_scores.begin(), sorted_scores.end(),
        [](const auto& a, const auto& b) {
            return a.second > b.second;
        }
    );
    
    // Return top N
    if (sorted_scores.size() > n) {
        sorted_scores.resize(n);
    }
    
    return sorted_scores;
}

} // namespace motif
