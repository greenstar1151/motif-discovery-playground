/**
 * @file kmer_counting.hpp
 * @brief K-mer Counting Operation
 * 
 * Core Operation 1: Extract all k-mers from sequences and count occurrences.
 * Uses ZOOPS (Zero Or One Per Sequence) model - each k-mer counted at most
 * once per sequence.
 * 
 * Complexity: O(N × L) where N = number of sequences, L = average length
 */

#ifndef MOTIF_KMER_COUNTING_HPP
#define MOTIF_KMER_COUNTING_HPP

#include "motif/types.hpp"

namespace motif {

/**
 * @brief Count k-mer occurrences using ZOOPS model
 * 
 * For each k-mer, counts the number of sequences it appears in
 * (not total occurrences). This is the standard model for motif discovery.
 * 
 * @param sequences List of DNA sequences
 * @param k K-mer length
 * @return Map of k-mer -> count (number of sequences containing it)
 * 
 * @example
 *   sequences = {"ACGTATATA", "TATATAGG", "AAACCCGGG"}
 *   k = 6
 *   result["TATATA"] = 2  // appears in sequences 0 and 1
 */
KmerCounts count_kmers(const SequenceList& sequences, size_t k);

/**
 * @brief Count k-mer occurrences (total count, not ZOOPS)
 * 
 * Alternative counting mode: counts every occurrence of each k-mer
 * across all sequences.
 * 
 * @param sequences List of DNA sequences
 * @param k K-mer length
 * @return Map of k-mer -> total occurrence count
 */
KmerCounts count_kmers_total(const SequenceList& sequences, size_t k);

} // namespace motif

#endif // MOTIF_KMER_COUNTING_HPP
