/**
 * @file fasta_reader.hpp
 * @brief FASTA File Reader for Real DNA Sequence Data
 * 
 * Utilities for loading sequences from FASTA format files.
 * Supports standard FASTA with multi-line sequences.
 */

#ifndef MOTIF_FASTA_READER_HPP
#define MOTIF_FASTA_READER_HPP

#include "motif/types.hpp"
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace motif {

/**
 * @brief Result of reading a FASTA file
 */
struct FastaData {
    SequenceList sequences;              ///< DNA sequences
    std::vector<std::string> headers;    ///< Sequence headers (names)
    
    size_t total_bases() const {
        size_t total = 0;
        for (const auto& seq : sequences) {
            total += seq.size();
        }
        return total;
    }
    
    double avg_length() const {
        if (sequences.empty()) return 0.0;
        return static_cast<double>(total_bases()) / sequences.size();
    }
};

/**
 * @brief Read sequences from a FASTA file
 * 
 * @param filepath Path to the FASTA file (.fa, .fasta)
 * @param max_sequences Maximum number of sequences to read (0 = unlimited)
 * @return FastaData containing sequences and headers
 * @throws std::runtime_error if file cannot be opened
 * 
 * @example
 *   auto data = read_fasta("upstream5000.fa");
 *   std::cout << "Loaded " << data.sequences.size() << " sequences\n";
 */
inline FastaData read_fasta(const std::string& filepath, size_t max_sequences = 0) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + filepath);
    }
    
    FastaData result;
    std::string line;
    std::string current_header;
    std::string current_seq;
    
    auto flush_sequence = [&]() {
        if (!current_seq.empty()) {
            // Normalize: uppercase, replace unknown bases with N
            for (char& c : current_seq) {
                c = std::toupper(static_cast<unsigned char>(c));
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                    c = 'N';
                }
            }
            result.sequences.push_back(std::move(current_seq));
            result.headers.push_back(std::move(current_header));
            current_seq.clear();
            current_header.clear();
        }
    };
    
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        
        // Remove trailing whitespace/carriage return
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        
        if (line[0] == '>') {
            // New sequence header
            flush_sequence();
            
            // Check limit
            if (max_sequences > 0 && result.sequences.size() >= max_sequences) {
                break;
            }
            
            // Parse header (everything after '>')
            current_header = line.substr(1);
        } else {
            // Sequence line - append to current
            current_seq += line;
        }
    }
    
    // Don't forget the last sequence
    if (max_sequences == 0 || result.sequences.size() < max_sequences) {
        flush_sequence();
    }
    
    return result;
}

/**
 * @brief Chop sequences into fixed-length segments
 * 
 * This is useful for preparing data for motif discovery,
 * especially when working with long sequences that need to be analyzed in smaller windows.
 * 
 * @param sequences Input sequences
 * @param segment_length Length of each segment
 * @param stride Step size between segments (default = segment_length for non-overlapping)
 * @return Chopped sequence list
 * 
 * @example
 *   // 5000bp -> 5 × 1000bp
 *   auto chopped = chop_sequences(data.sequences, 1000);
 */
inline SequenceList chop_sequences(
    const SequenceList& sequences,
    size_t segment_length,
    size_t stride = 0
) {
    if (stride == 0) stride = segment_length;
    
    SequenceList result;
    result.reserve(sequences.size() * 2);  // Rough estimate
    
    for (const auto& seq : sequences) {
        if (seq.length() < segment_length) {
            // Keep short sequences as-is
            result.push_back(seq);
            continue;
        }
        
        for (size_t i = 0; i + segment_length <= seq.length(); i += stride) {
            result.push_back(seq.substr(i, segment_length));
        }
    }
    
    return result;
}

/**
 * @brief Sample a subset of sequences randomly
 * 
 * @param sequences Input sequences
 * @param count Number of sequences to sample
 * @param seed Random seed
 * @return Sampled sequence list
 */
inline SequenceList sample_sequences(
    const SequenceList& sequences,
    size_t count,
    unsigned int seed = 42
) {
    if (count >= sequences.size()) {
        return sequences;
    }
    
    std::mt19937 rng(seed);
    std::vector<size_t> indices(sequences.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    
    SequenceList result;
    result.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        result.push_back(sequences[indices[i]]);
    }
    
    return result;
}

/**
 * @brief Compute base composition statistics
 */
struct BaseComposition {
    size_t count_A = 0;
    size_t count_C = 0;
    size_t count_G = 0;
    size_t count_T = 0;
    size_t count_N = 0;
    
    size_t total() const { return count_A + count_C + count_G + count_T + count_N; }
    double gc_content() const {
        size_t acgt = count_A + count_C + count_G + count_T;
        return acgt > 0 ? static_cast<double>(count_C + count_G) / acgt : 0.0;
    }
    
    double freq_A() const { return total() > 0 ? static_cast<double>(count_A) / total() : 0.0; }
    double freq_C() const { return total() > 0 ? static_cast<double>(count_C) / total() : 0.0; }
    double freq_G() const { return total() > 0 ? static_cast<double>(count_G) / total() : 0.0; }
    double freq_T() const { return total() > 0 ? static_cast<double>(count_T) / total() : 0.0; }
};

inline BaseComposition compute_base_composition(const SequenceList& sequences) {
    BaseComposition comp;
    for (const auto& seq : sequences) {
        for (char c : seq) {
            switch (std::toupper(static_cast<unsigned char>(c))) {
                case 'A': ++comp.count_A; break;
                case 'C': ++comp.count_C; break;
                case 'G': ++comp.count_G; break;
                case 'T': ++comp.count_T; break;
                default:  ++comp.count_N; break;
            }
        }
    }
    return comp;
}

// =============================================================================
// Binding Site Injection (for real data benchmarking)
// =============================================================================

/**
 * @brief Inject real binding sites into baseline sequences
 * 
 * Take actual binding site sequences from a .sites file (FASTA format)
 * and insert them at random positions into baseline sequences.
 * 
 * @param baseline_sequences Background sequences (e.g., from upstream5000.fa)
 * @param binding_sites Real binding site sequences (e.g., from MA0007.2.sites)
 * @param seed Random seed for reproducibility
 * @return Sequences with one binding site inserted per sequence
 * 
 * @example
 *   auto baseline = read_fasta("upstream5000.fa");
 *   auto sites = read_fasta("MA0007.2.sites");
 *   auto chopped = chop_sequences(baseline.sequences, 1000);
 *   auto sampled_baseline = sample_sequences(chopped, 32768);
 *   auto sampled_sites = sample_sequences(sites.sequences, 32768);
 *   auto primary = inject_binding_sites(sampled_baseline, sampled_sites);
 */
inline SequenceList inject_binding_sites(
    const SequenceList& baseline_sequences,
    const SequenceList& binding_sites,
    unsigned int seed = 42
) {
    if (baseline_sequences.empty() || binding_sites.empty()) {
        return baseline_sequences;
    }
    
    std::mt19937 rng(seed);
    
    SequenceList result;
    result.reserve(baseline_sequences.size());
    
    for (size_t i = 0; i < baseline_sequences.size(); ++i) {
        const auto& baseline = baseline_sequences[i];
        // Cycle through binding sites if fewer than baseline sequences
        const auto& site = binding_sites[i % binding_sites.size()];
        
        if (baseline.length() < site.length()) {
            // Baseline too short, keep as-is
            result.push_back(baseline);
            continue;
        }
        
        // Random insertion position
        std::uniform_int_distribution<size_t> pos_dist(0, baseline.length() - site.length());
        size_t pos = pos_dist(rng);
        
        // Insert binding site
        std::string new_seq = baseline;
        new_seq.replace(pos, site.length(), site);
        result.push_back(std::move(new_seq));
    }
    
    return result;
}

/**
 * @brief Inject binding sites with configurable injection rate
 * 
 * @param baseline_sequences Background sequences
 * @param binding_sites Real binding site sequences
 * @param injection_rate Probability of injection per sequence (0.0-1.0)
 * @param seed Random seed
 * @return Sequences with binding sites injected at given rate
 */
inline SequenceList inject_binding_sites_with_rate(
    const SequenceList& baseline_sequences,
    const SequenceList& binding_sites,
    double injection_rate,
    unsigned int seed = 42
) {
    if (baseline_sequences.empty() || binding_sites.empty()) {
        return baseline_sequences;
    }
    
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    
    SequenceList result;
    result.reserve(baseline_sequences.size());
    
    size_t site_idx = 0;
    for (const auto& baseline : baseline_sequences) {
        if (prob_dist(rng) < injection_rate && !binding_sites.empty()) {
            const auto& site = binding_sites[site_idx % binding_sites.size()];
            ++site_idx;
            
            if (baseline.length() >= site.length()) {
                std::uniform_int_distribution<size_t> pos_dist(0, baseline.length() - site.length());
                size_t pos = pos_dist(rng);
                
                std::string new_seq = baseline;
                new_seq.replace(pos, site.length(), site);
                result.push_back(std::move(new_seq));
                continue;
            }
        }
        result.push_back(baseline);
    }
    
    return result;
}

/**
 * @brief Generate primary and control sets from real data
 * 
 * 1. Load baseline sequences
 * 2. Chop into 1000bp segments
 * 3. Sample N sequences for primary set
 * 4. Sample M binding sites and inject into primary set
 * 5. Sample N different sequences for control set (no injection)
 * 
 * @param baseline_file Path to baseline FASTA file
 * @param sites_file Path to binding sites FASTA file
 * @param num_sequences Number of sequences per set
 * @param segment_length Length to chop baseline sequences
 * @param seed Random seed
 * @return Pair of (primary_sequences, control_sequences)
 */
inline std::pair<SequenceList, SequenceList> generate_real_test_data(
    const std::string& baseline_file,
    const std::string& sites_file,
    size_t num_sequences,
    size_t segment_length = 1000,
    unsigned int seed = 42
) {
    // Load baseline sequences
    auto baseline_data = read_fasta(baseline_file);
    auto chopped = chop_sequences(baseline_data.sequences, segment_length);
    
    // Load binding sites
    auto sites_data = read_fasta(sites_file);
    
    // Sample for primary set
    auto primary_baseline = sample_sequences(chopped, num_sequences, seed);
    auto sampled_sites = sample_sequences(sites_data.sequences, num_sequences, seed + 1);
    auto primary = inject_binding_sites(primary_baseline, sampled_sites, seed + 2);
    
    // Sample different sequences for control set (use different seed region)
    auto control = sample_sequences(chopped, num_sequences, seed + 1000);
    
    return {primary, control};
}

} // namespace motif

#endif // MOTIF_FASTA_READER_HPP
