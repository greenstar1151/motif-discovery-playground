/**
 * @file kmer_counting.cpp
 * @brief Implementation of K-mer Counting (Core Operation 1)
 */

#include "motif/kmer_counting.hpp"
#include <unordered_set>

namespace motif {

KmerCounts count_kmers(const SequenceList& sequences, size_t k) {
    KmerCounts counts;
    
    // Process each sequence
    for (const auto& seq : sequences) {
        // Skip sequences shorter than k
        if (seq.length() < k) {
            continue;
        }
        
        // Track k-mers seen in this sequence (ZOOPS model)
        std::unordered_set<std::string> seen_in_seq;
        seen_in_seq.reserve(seq.length() - k + 1);
        
        // Slide window over sequence
        for (size_t i = 0; i <= seq.length() - k; ++i) {
            std::string kmer = seq.substr(i, k);
            
            // Skip k-mers containing masked bases
            if (!is_valid_kmer(kmer)) {
                continue;
            }
            
            // ZOOPS: count each k-mer only once per sequence
            if (seen_in_seq.find(kmer) == seen_in_seq.end()) {
                counts[kmer]++;
                seen_in_seq.insert(kmer);
            }
        }
    }
    
    return counts;
}

KmerCounts count_kmers_total(const SequenceList& sequences, size_t k) {
    KmerCounts counts;
    
    // Process each sequence
    for (const auto& seq : sequences) {
        // Skip sequences shorter than k
        if (seq.length() < k) {
            continue;
        }
        
        // Slide window over sequence
        for (size_t i = 0; i <= seq.length() - k; ++i) {
            std::string kmer = seq.substr(i, k);
            
            // Skip k-mers containing masked bases
            if (!is_valid_kmer(kmer)) {
                continue;
            }
            
            // Count every occurrence
            counts[kmer]++;
        }
    }
    
    return counts;
}

} // namespace motif
