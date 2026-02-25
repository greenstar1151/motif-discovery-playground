/**
 * @file kmer_counting.cpp
 * @brief Implementation of K-mer counting primitives.
 */

#include "motif/kmer_counting.hpp"
#include <unordered_set>

namespace motif {

KmerCounts count_kmers(const SequenceList& sequences, size_t k) {
    KmerCounts counts;
    if (k == 0) {
        return counts;
    }

    for (const auto& seq : sequences) {
        if (seq.size() < k) {
            continue;
        }

        std::unordered_set<std::string> seen_in_seq;
        seen_in_seq.reserve(seq.size() - k + 1);

        for (size_t pos = 0; pos + k <= seq.size(); ++pos) {
            const std::string kmer = seq.substr(pos, k);
            if (!is_valid_kmer(kmer)) {
                continue;
            }

            if (seen_in_seq.insert(kmer).second) {
                counts[kmer] += 1;
            }
        }
    }

    return counts;
}

KmerCounts count_kmers_total(const SequenceList& sequences, size_t k) {
    KmerCounts counts;
    if (k == 0) {
        return counts;
    }

    for (const auto& seq : sequences) {
        if (seq.size() < k) {
            continue;
        }

        for (size_t pos = 0; pos + k <= seq.size(); ++pos) {
            const std::string kmer = seq.substr(pos, k);
            if (!is_valid_kmer(kmer)) {
                continue;
            }
            counts[kmer] += 1;
        }
    }

    return counts;
}

} // namespace motif
