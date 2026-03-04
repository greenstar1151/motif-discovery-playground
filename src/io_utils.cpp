/**
 * @file io_utils.cpp
 * @brief Implementation of I/O and data-preparation utilities.
 */

#include "motif/io_utils.hpp"
#include "motif/markov.hpp"  // kMinProb

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace motif {

// ─────────────────────────────────────────────────────────────────────────────
// shuffle_sequence_kmer_preserving
// ─────────────────────────────────────────────────────────────────────────────

std::string shuffle_sequence_kmer_preserving(
    const std::string& seq,
    std::mt19937& rng,
    size_t k
) {
    if (seq.size() <= k) {
        return seq;
    }

    if (k <= 1) {
        std::string shuffled = seq;
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        return shuffled;
    }

    const size_t prefix_len = k - 1;

    // Build the de-Bruijn multi-graph
    std::map<std::string, std::vector<std::pair<char, size_t>>> adj;
    std::vector<bool> used;
    size_t edge_count = 0;

    for (size_t i = 0; i + k <= seq.size(); ++i) {
        const std::string vertex = seq.substr(i, prefix_len);
        const char next_char = seq[i + prefix_len];
        adj[vertex].push_back({next_char, edge_count++});
        used.push_back(false);
    }

    // Random-shuffle adjacency lists
    for (auto& [_, neighbors] : adj) {
        std::shuffle(neighbors.begin(), neighbors.end(), rng);
    }

    // Hierholzer Euler-tour
    const std::string start_vertex = seq.substr(0, prefix_len);
    std::vector<char> path_chars;
    std::vector<std::string> stack;
    stack.push_back(start_vertex);

    std::map<std::string, size_t> adj_pos;
    for (const auto& [v, _] : adj) {
        adj_pos[v] = 0;
    }

    while (!stack.empty()) {
        const std::string v = stack.back();

        bool found = false;
        while (adj_pos.count(v) && adj_pos[v] < adj[v].size()) {
            auto [next_char, edge_idx] = adj[v][adj_pos[v]];
            adj_pos[v] += 1;

            if (!used[edge_idx]) {
                used[edge_idx] = true;
                const std::string next_vertex = v.substr(1) + next_char;
                stack.push_back(next_vertex);
                found = true;
                break;
            }
        }

        if (!found) {
            if (stack.size() > 1) {
                path_chars.push_back(v.back());
            }
            stack.pop_back();
        }
    }

    std::string result = start_vertex;
    for (auto it = path_chars.rbegin(); it != path_chars.rend(); ++it) {
        result.push_back(*it);
    }

    // Fallback: if Euler tour did not cover all edges, simple shuffle
    if (result.size() != seq.size()) {
        result = seq;
        std::shuffle(result.begin(), result.end(), rng);
    }

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// generate_control_sequences
// ─────────────────────────────────────────────────────────────────────────────

SequenceList generate_control_sequences(
    const SequenceList& input,
    unsigned int seed,
    size_t k
) {
    std::mt19937 rng(seed);
    SequenceList control;
    control.reserve(input.size());

    for (const auto& seq : input) {
        control.push_back(shuffle_sequence_kmer_preserving(seq, rng, k));
    }

    return control;
}

// ─────────────────────────────────────────────────────────────────────────────
// load_jaspar_pwm
// ─────────────────────────────────────────────────────────────────────────────

PWM load_jaspar_pwm(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open JASPAR file: " + filepath);
    }

    PWM pwm;
    std::string line;
    std::array<std::vector<double>, 4> counts;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') {
            continue;
        }

        const size_t base_idx = base_to_index(line[0]);
        if (base_idx == SIZE_MAX) {
            continue;
        }

        const size_t start = line.find('[');
        const size_t end   = line.find(']');
        if (start == std::string::npos || end == std::string::npos ||
            end <= start + 1) {
            continue;
        }

        std::istringstream iss(line.substr(start + 1, end - start - 1));
        double value = 0.0;
        while (iss >> value) {
            counts[base_idx].push_back(value);
        }
    }

    if (counts[0].empty()) {
        throw std::runtime_error("Invalid or empty JASPAR file");
    }

    pwm.width = counts[0].size();
    if (pwm.width > MAX_MOTIF_WIDTH) {
        throw std::runtime_error("Motif width exceeds MAX_MOTIF_WIDTH");
    }

    for (size_t pos = 0; pos < pwm.width; ++pos) {
        double total = 0.0;
        for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
            total += counts[b][pos];
        }

        const double pseudo = 0.1;
        for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
            pwm.matrix[b][pos] =
                (counts[b][pos] + pseudo) /
                std::max(total + 4.0 * pseudo, kMinProb);
        }
    }

    return pwm;
}

} // namespace motif
