/**
 * @file morbius_eval.cpp
 * @brief Morbius Benchmark Evaluation - Generate result TSV for accuracy comparison
 * 
 * This tool processes a FASTA dataset and outputs motif predictions in TSV format
 * compatible with the Morbius benchmark evaluation protocol.
 * 
 * Output format (TSV):
 *   seq_id    position    motif
 *   0         859         TATATA
 *   1         729         ATATAT
 *   ...
 * 
 * Usage:
 *   ./morbius_eval <input.fasta> <output.tsv> [options]
 * 
 * Options:
 *   -k <int>              K-mer length for seed discovery (default: 16)
 *   -m <int>              Motif width for PWM scoring (default: same as k)
 *   --control <file>      Control FASTA file (default: use shuffled input)
 *   --seed <string>       Use specific seed motif instead of discovering
 *   --pwm <file>          Load PWM from JASPAR file
 *   --top <int>           Number of top seeds to try (default: 1)
 *   --help                Show this help
 */

#include "motif.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <random>

using namespace motif;
using namespace std::chrono;

// =============================================================================
// Utilities
// =============================================================================

void print_usage(const char* prog) {
    std::cerr << "\nUsage: " << prog << " <input.fasta> <output.tsv> [options]\n\n"
              << "Options:\n"
              << "  -k <int>              K-mer length for seed discovery (default: 16)\n"
              << "  -m <int>              Motif width for PWM scoring (default: same as k)\n"
              << "  --control <file>      Control FASTA file (default: use shuffled input)\n"
              << "  --seed <string>       Use specific seed motif instead of discovering\n"
              << "  --top <int>           Number of top seeds to try (default: 1)\n"
              << "  --threshold <float>   PWM score threshold (default: 0.0, use best match)\n"
              << "  --help                Show this help\n\n"
              << "Output format (TSV):\n"
              << "  seq_id    position    motif\n\n"
              << "Example:\n"
              << "  " << prog << " DATASET_DNA_3.fasta DATASET_DNA_3_result.tsv -k 16\n";
}

/**
 * @brief Shuffle a sequence preserving k-mer (default: dinucleotide) frequencies
 * 
 * Uses the Euler path method with Hierholzer's algorithm.
 * For k=2 (dinucleotide), this preserves:
 * - Single nucleotide frequencies (A, C, G, T counts)
 * - Dinucleotide frequencies (AA, AC, AG, ..., TT counts)
 * - First and last nucleotides of the sequence
 * 
 * Based on Altschul & Erickson (1985) and ushuffle algorithm.
 */
std::string shuffle_sequence_kmer_preserving(const std::string& seq, std::mt19937& rng, size_t k = 2) {
    if (seq.length() <= k) {
        return seq;
    }
    
    // For k=1, just do simple shuffle
    if (k == 1) {
        std::string shuffled = seq;
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        return shuffled;
    }
    
    // For dinucleotide (k=2), use graph-based Euler path
    // Vertices = (k-1)-mers, Edges = k-mers
    const size_t prefix_len = k - 1;
    
    // Build adjacency list: vertex -> list of (next_char, edge_index)
    std::map<std::string, std::vector<std::pair<char, size_t>>> adj;
    std::vector<bool> used;
    size_t edge_count = 0;
    
    for (size_t i = 0; i + k <= seq.length(); ++i) {
        std::string vertex = seq.substr(i, prefix_len);
        char next_char = seq[i + prefix_len];
        adj[vertex].push_back({next_char, edge_count++});
        used.push_back(false);
    }
    
    // Shuffle adjacency lists for randomization
    for (auto& [vertex, neighbors] : adj) {
        std::shuffle(neighbors.begin(), neighbors.end(), rng);
    }
    
    // Hierholzer's algorithm to find Euler path
    std::string start_vertex = seq.substr(0, prefix_len);
    std::vector<char> path_chars;
    std::vector<std::string> stack;
    stack.push_back(start_vertex);
    
    // Track current position in each adjacency list
    std::map<std::string, size_t> adj_pos;
    for (auto& [v, _] : adj) {
        adj_pos[v] = 0;
    }
    
    while (!stack.empty()) {
        std::string v = stack.back();
        
        // Find next unused edge from v
        bool found = false;
        while (adj_pos.count(v) && adj_pos[v] < adj[v].size()) {
            auto& [next_char, edge_idx] = adj[v][adj_pos[v]];
            adj_pos[v]++;
            
            if (!used[edge_idx]) {
                used[edge_idx] = true;
                std::string next_vertex = v.substr(1) + next_char;
                stack.push_back(next_vertex);
                found = true;
                break;
            }
        }
        
        if (!found) {
            // No more edges from v, add to path
            if (stack.size() > 1) {
                // Extract the character that led to this vertex
                std::string prev = stack[stack.size() - 2];
                // The edge from prev to v determines the character
                // v = prev[1:] + char, so char = v.back()
                path_chars.push_back(v.back());
            }
            stack.pop_back();
        }
    }
    
    // Build result: start with first (k-1) chars, then add path in reverse
    std::string result = start_vertex;
    for (auto it = path_chars.rbegin(); it != path_chars.rend(); ++it) {
        result += *it;
    }
    
    // Verify we used all edges
    if (result.length() != seq.length()) {
        // Fallback to simple shuffle if Euler path failed
        // (can happen with disconnected graph or non-Eulerian)
        result = seq;
        std::shuffle(result.begin(), result.end(), rng);
    }
    
    return result;
}

/**
 * @brief Simple shuffle (does not preserve k-mer frequencies)
 */
std::string shuffle_sequence_simple(const std::string& seq, std::mt19937& rng) {
    std::string shuffled = seq;
    std::shuffle(shuffled.begin(), shuffled.end(), rng);
    return shuffled;
}

/**
 * @brief Generate control sequences by shuffling input sequences
 * Uses k-mer preserving shuffle (default k=2 for dinucleotide preservation)
 */
SequenceList generate_control_sequences(const SequenceList& input, unsigned int seed, size_t k = 2) {
    std::mt19937 rng(seed);
    SequenceList control;
    control.reserve(input.size());
    
    for (const auto& seq : input) {
        control.push_back(shuffle_sequence_kmer_preserving(seq, rng, k));
    }
    
    return control;
}

/**
 * @brief Parse JASPAR PWM file
 * Format:
 *   >MA0007.2 AR
 *   A [ count1 count2 ... ]
 *   C [ count1 count2 ... ]
 *   G [ count1 count2 ... ]
 *   T [ count1 count2 ... ]
 */
PWM load_jaspar_pwm(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JASPAR file: " + filepath);
    }
    
    PWM pwm;
    std::string line;
    std::array<std::vector<double>, 4> counts;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        
        // Parse base line: A [ 1 2 3 4 ... ]
        char base = line[0];
        size_t base_idx = base_to_index(base);
        if (base_idx == SIZE_MAX) continue;
        
        // Find values between [ ]
        size_t start = line.find('[');
        size_t end = line.find(']');
        if (start == std::string::npos || end == std::string::npos) continue;
        
        std::string values = line.substr(start + 1, end - start - 1);
        std::istringstream iss(values);
        double val;
        while (iss >> val) {
            counts[base_idx].push_back(val);
        }
    }
    
    // Convert counts to probabilities
    if (counts[0].empty()) {
        throw std::runtime_error("Empty or invalid JASPAR file");
    }
    
    pwm.width = counts[0].size();
    if (pwm.width > MAX_MOTIF_WIDTH) {
        throw std::runtime_error("Motif width exceeds maximum: " + std::to_string(pwm.width));
    }
    
    for (size_t pos = 0; pos < pwm.width; ++pos) {
        double total = 0;
        for (size_t b = 0; b < 4; ++b) {
            total += counts[b][pos];
        }
        
        // Add pseudocount and normalize
        double pseudo = 0.1;
        for (size_t b = 0; b < 4; ++b) {
            pwm.set(b, pos, (counts[b][pos] + pseudo) / (total + 4 * pseudo));
        }
    }
    
    return pwm;
}

// =============================================================================
// Main Processing
// =============================================================================

struct Config {
    std::string input_file;
    std::string output_file;
    std::string control_file;
    std::string seed_motif;
    std::string jaspar_file;
    size_t k = 16;
    size_t motif_width = 0;  // 0 means use k
    size_t top_seeds = 1;
    double threshold = 0.0;  // 0 means use best match per sequence
    unsigned int rng_seed = 42;
};

struct MotifHit {
    size_t seq_id;
    int position;
    std::string motif;
    double score;
};

/**
 * @brief Run motif discovery and generate results
 */
std::vector<MotifHit> run_discovery(
    const FastaData& input_data,
    const SequenceList& control_seqs,
    const Config& config
) {
    const size_t k = config.k;
    (void)config.motif_width;  // Reserved for future use
    
    std::vector<MotifHit> results;
    results.reserve(input_data.sequences.size());
    
    // =========================================================================
    // Step 1: Discover or use provided seed
    // =========================================================================
    PWM pwm;
    std::string best_seed;
    
    if (!config.jaspar_file.empty()) {
        // Load PWM from JASPAR file
        std::cerr << "[Info] Loading PWM from JASPAR file: " << config.jaspar_file << std::endl;
        pwm = load_jaspar_pwm(config.jaspar_file);
        best_seed = get_consensus(pwm);
        std::cerr << "[Info] PWM consensus: " << best_seed << " (width: " << pwm.width << ")" << std::endl;
    } else if (!config.seed_motif.empty()) {
        // Use provided seed
        best_seed = config.seed_motif;
        std::cerr << "[Info] Using provided seed: " << best_seed << std::endl;
        pwm = build_pwm_from_seed(best_seed, input_data.sequences);
    } else {
        // Discover seed from k-mer enrichment
        std::cerr << "[Info] Discovering seed with k=" << k << std::endl;
        
        auto pos_counts = count_kmers(input_data.sequences, k);
        auto neg_counts = count_kmers(control_seqs, k);
        
        auto scores = calculate_enrichment_scores(
            pos_counts, neg_counts,
            static_cast<int>(input_data.sequences.size()),
            static_cast<int>(control_seqs.size())
        );
        
        auto top_seeds = get_top_seeds(scores, config.top_seeds);
        
        if (top_seeds.empty()) {
            std::cerr << "[Error] No k-mers found!" << std::endl;
            return results;
        }
        
        best_seed = top_seeds[0].first;
        double best_score = top_seeds[0].second;
        
        std::cerr << "[Info] Top seed: " << best_seed 
                  << " (enrichment: " << std::fixed << std::setprecision(4) << best_score << ")" << std::endl;
        
        if (config.top_seeds > 1 && top_seeds.size() > 1) {
            std::cerr << "[Info] Other top seeds:" << std::endl;
            for (size_t i = 1; i < std::min(config.top_seeds, top_seeds.size()); ++i) {
                std::cerr << "       " << top_seeds[i].first 
                          << " (" << std::fixed << std::setprecision(4) << top_seeds[i].second << ")" << std::endl;
            }
        }
        
        // Build PWM from seed
        pwm = build_pwm_from_seed(best_seed, input_data.sequences);
    }
    
    std::cerr << "[Info] PWM consensus: " << get_consensus(pwm) << std::endl;
    
    // =========================================================================
    // Step 2: Scan each sequence and find best match
    // =========================================================================
    std::cerr << "[Info] Scanning " << input_data.sequences.size() << " sequences..." << std::endl;
    
    for (size_t i = 0; i < input_data.sequences.size(); ++i) {
        const auto& seq = input_data.sequences[i];
        
        // Find best match position using PWM scoring
        auto [best_pos, best_score] = find_best_match(seq, pwm);
        
        MotifHit hit;
        hit.seq_id = i;
        hit.score = best_score;
        
        if (best_pos >= 0 && (config.threshold == 0.0 || best_score >= config.threshold)) {
            hit.position = best_pos;
            // Extract the motif at the predicted position
            if (static_cast<size_t>(best_pos) + pwm.width <= seq.length()) {
                hit.motif = seq.substr(best_pos, pwm.width);
            } else {
                hit.motif = seq.substr(best_pos);
            }
        } else {
            // No significant match found
            hit.position = -1;
            hit.motif = "";
        }
        
        results.push_back(hit);
        
        // Progress indicator
        if ((i + 1) % 10000 == 0) {
            std::cerr << "[Progress] " << (i + 1) << " / " << input_data.sequences.size() << std::endl;
        }
    }
    
    return results;
}

/**
 * @brief Write results to TSV file
 */
void write_results_tsv(
    const std::string& filepath,
    const std::vector<MotifHit>& results,
    const std::vector<std::string>& headers
) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filepath);
    }
    
    // Write header
    out << "seq_id\tposition\tmotif\n";
    
    // Write results
    for (const auto& hit : results) {
        // Try to parse seq_id from header if available, otherwise use index
        std::string seq_id_str;
        if (hit.seq_id < headers.size()) {
            // Use header as seq_id
            // For Morbius format, header is just a number like "0", "1", etc.
            // For other formats, extract the first word (before any space/tab)
            const std::string& header = headers[hit.seq_id];
            size_t space_pos = header.find_first_of(" \t");
            if (space_pos != std::string::npos) {
                seq_id_str = header.substr(0, space_pos);
            } else {
                seq_id_str = header;
            }
        } else {
            seq_id_str = std::to_string(hit.seq_id);
        }
        
        if (hit.position >= 0) {
            out << seq_id_str << "\t" << hit.position << "\t" << hit.motif << "\n";
        } else {
            // No match found - output with -1 position and empty motif
            out << seq_id_str << "\t-1\t\n";
        }
    }
    
    out.close();
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }
    
    Config config;
    
    // Parse arguments
    int i = 1;
    while (i < argc) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-k" && i + 1 < argc) {
            config.k = std::stoull(argv[++i]);
        } else if (arg == "-m" && i + 1 < argc) {
            config.motif_width = std::stoull(argv[++i]);
        } else if (arg == "--control" && i + 1 < argc) {
            config.control_file = argv[++i];
        } else if (arg == "--seed" && i + 1 < argc) {
            config.seed_motif = argv[++i];
        } else if (arg == "--pwm" && i + 1 < argc) {
            config.jaspar_file = argv[++i];
        } else if (arg == "--top" && i + 1 < argc) {
            config.top_seeds = std::stoull(argv[++i]);
        } else if (arg == "--threshold" && i + 1 < argc) {
            config.threshold = std::stod(argv[++i]);
        } else if (arg == "--rng-seed" && i + 1 < argc) {
            config.rng_seed = static_cast<unsigned int>(std::stoi(argv[++i]));
        } else if (config.input_file.empty()) {
            config.input_file = arg;
        } else if (config.output_file.empty()) {
            config.output_file = arg;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
        ++i;
    }
    
    if (config.input_file.empty() || config.output_file.empty()) {
        std::cerr << "Error: Input and output files are required.\n";
        print_usage(argv[0]);
        return 1;
    }
    
    try {
        auto start_time = high_resolution_clock::now();
        
        // =====================================================================
        // Load input data
        // =====================================================================
        std::cerr << "=========================================================\n";
        std::cerr << "  Morbius Benchmark Evaluation - Motif Discovery\n";
        std::cerr << "=========================================================\n\n";
        
        std::cerr << "[Info] Loading input: " << config.input_file << std::endl;
        auto input_data = read_fasta(config.input_file);
        std::cerr << "[Info] Loaded " << input_data.sequences.size() << " sequences" << std::endl;
        
        if (input_data.sequences.empty()) {
            std::cerr << "[Error] No sequences found in input file!" << std::endl;
            return 1;
        }
        
        std::cerr << "[Info] Average sequence length: " << std::fixed << std::setprecision(1) 
                  << input_data.avg_length() << " bp" << std::endl;
        
        // =====================================================================
        // Prepare control sequences
        // =====================================================================
        SequenceList control_seqs;
        
        if (!config.control_file.empty()) {
            std::cerr << "[Info] Loading control: " << config.control_file << std::endl;
            auto control_data = read_fasta(config.control_file);
            control_seqs = std::move(control_data.sequences);
            std::cerr << "[Info] Loaded " << control_seqs.size() << " control sequences" << std::endl;
        } else {
            std::cerr << "[Info] Generating control sequences by shuffling..." << std::endl;
            control_seqs = generate_control_sequences(input_data.sequences, config.rng_seed);
        }
        
        // =====================================================================
        // Run discovery
        // =====================================================================
        std::cerr << std::endl;
        auto results = run_discovery(input_data, control_seqs, config);
        
        // =====================================================================
        // Write results
        // =====================================================================
        std::cerr << std::endl;
        std::cerr << "[Info] Writing results to: " << config.output_file << std::endl;
        write_results_tsv(config.output_file, results, input_data.headers);
        
        // =====================================================================
        // Summary
        // =====================================================================
        auto end_time = high_resolution_clock::now();
        double elapsed = duration_cast<milliseconds>(end_time - start_time).count() / 1000.0;
        
        int hits_found = 0;
        for (const auto& hit : results) {
            if (hit.position >= 0) ++hits_found;
        }
        
        std::cerr << std::endl;
        std::cerr << "=========================================================\n";
        std::cerr << "  Summary\n";
        std::cerr << "=========================================================\n";
        std::cerr << "  Total sequences : " << results.size() << std::endl;
        std::cerr << "  Motifs found    : " << hits_found << std::endl;
        std::cerr << "  Coverage        : " << std::fixed << std::setprecision(1) 
                  << (100.0 * hits_found / results.size()) << "%" << std::endl;
        std::cerr << "  Elapsed time    : " << std::fixed << std::setprecision(2) 
                  << elapsed << " seconds" << std::endl;
        std::cerr << "  Output file     : " << config.output_file << std::endl;
        std::cerr << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "[Error] " << e.what() << std::endl;
        return 1;
    }
}
