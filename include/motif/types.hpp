/**
 * @file types.hpp
 * @brief Common type definitions for motif discovery
 * 
 * Defines the fundamental data structures used across all operations.
 * Designed for clarity and ease of hardware mapping.
 */

#ifndef MOTIF_TYPES_HPP
#define MOTIF_TYPES_HPP

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace motif {

// =============================================================================
// Constants
// =============================================================================

/// DNA alphabet size (A, C, G, T)
constexpr size_t ALPHABET_SIZE = 4;

/// Maximum supported motif width
constexpr size_t MAX_MOTIF_WIDTH = 32;

/// Base indices for consistent mapping
constexpr size_t BASE_A = 0;
constexpr size_t BASE_C = 1;
constexpr size_t BASE_G = 2;
constexpr size_t BASE_T = 3;

/// Masking character
constexpr char MASK_CHAR = 'N';

// =============================================================================
// Type Aliases
// =============================================================================

/// A DNA sequence represented as a string
using Sequence = std::string;

/// Collection of DNA sequences
using SequenceList = std::vector<Sequence>;

/// K-mer count map: kmer string -> count
using KmerCounts = std::map<std::string, int>;

/// Enrichment score map: kmer string -> score
using EnrichmentScores = std::map<std::string, double>;

// =============================================================================
// Position Weight Matrix (PWM)
// =============================================================================

/**
 * @struct PWM
 * @brief Position Weight Matrix for motif representation
 * 
 * A 4 x W matrix where each row represents a nucleotide (A, C, G, T)
 * and each column represents a position in the motif.
 * 
 * Values are probabilities (0.0 to 1.0) that sum to 1.0 at each position.
 */
struct PWM {
    size_t width = 0;                                    ///< Motif width (number of positions)
    std::array<std::array<double, MAX_MOTIF_WIDTH>, 
               ALPHABET_SIZE> matrix = {};               ///< [base][position] -> probability
    
    /// Get probability for a specific base at a position
    double at(size_t base_idx, size_t pos) const {
        return matrix[base_idx][pos];
    }
    
    /// Set probability for a specific base at a position
    void set(size_t base_idx, size_t pos, double prob) {
        matrix[base_idx][pos] = prob;
    }
    
    /// Reset the PWM
    void clear() {
        width = 0;
        for (auto& row : matrix) {
            row.fill(0.0);
        }
    }
};

// =============================================================================
// Discovered Motif Result
// =============================================================================

/**
 * @struct DiscoveredMotif
 * @brief Result of motif discovery containing all relevant information
 */
struct DiscoveredMotif {
    std::string consensus;      ///< Best k-mer (consensus sequence)
    PWM pwm;                    ///< Position Weight Matrix
    double score = 0.0;         ///< Enrichment score
    int count_primary = 0;      ///< Count in primary set
    int count_control = 0;      ///< Count in control set
};

// =============================================================================
// Utility Functions
// =============================================================================

/**
 * @brief Convert nucleotide character to index
 * @param base Nucleotide character (A, C, G, T)
 * @return Index (0-3) or SIZE_MAX if invalid
 */
inline size_t base_to_index(char base) {
    switch (base) {
        case 'A': case 'a': return BASE_A;
        case 'C': case 'c': return BASE_C;
        case 'G': case 'g': return BASE_G;
        case 'T': case 't': return BASE_T;
        default: return SIZE_MAX;  // Invalid base
    }
}

/**
 * @brief Convert index to nucleotide character
 * @param idx Index (0-3)
 * @return Nucleotide character
 */
inline char index_to_base(size_t idx) {
    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    return (idx < ALPHABET_SIZE) ? bases[idx] : 'N';
}

/**
 * @brief Check if a k-mer contains only valid bases (no N)
 * @param kmer K-mer string
 * @return true if all bases are valid (A, C, G, T)
 */
inline bool is_valid_kmer(const std::string& kmer) {
    for (char c : kmer) {
        if (base_to_index(c) == SIZE_MAX) {
            return false;
        }
    }
    return true;
}

} // namespace motif

#endif // MOTIF_TYPES_HPP
