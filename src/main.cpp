/**
 * @file main.cpp
 * @brief Integrated software reference pipeline for motif discovery.
 */

#include "motif.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>

using namespace motif;

namespace {

constexpr double kMinProb = 1e-12;

struct PipelineConfig {
    bool use_real_data = false;
    std::string baseline_file = "data/upstream5000.fa";
    std::string sites_file = "data/MA0007.2.sites";

    size_t num_sequences = 12000;
    size_t seq_length = 100;
    size_t segment_length = 1000;
    std::string synthetic_motif = "TATATA";
    double synthetic_injection_rate = 0.45;
    unsigned int seed = 42;

    std::vector<size_t> k_values = {6, 8, 9, 10, 12};
    double significant_z = 6.0;
    double sub_significant_z = 3.5;
    size_t top_seeds_per_k = 25;
    size_t max_seed_trials_per_round = 80;

    int max_motifs = 5;
    int max_em_iterations = 20;
    int em_patience = 2;
    double em_min_improvement = 1e-3;
    size_t min_sites_for_em = 8;
    size_t max_sites_for_mstep = 60000;
    double pseudocount = 0.2;
    double min_enrichment_to_report = 1.20;
    double max_mask_fraction = 0.85;
};

struct SeedCandidate {
    std::string kmer;
    size_t width = 0;
    double z_score = 0.0;
    int pos_count = 0;
    int neg_count = 0;
};

struct MatchResult {
    int position = -1;
    double score = -std::numeric_limits<double>::infinity();
};

struct ThresholdResult {
    bool valid = false;
    double threshold = 0.0;
    double enrichment = 0.0;
    int primary_hits = 0;
    int control_hits = 0;
};

struct MotifRoundResult {
    std::string seed;
    std::string consensus;
    PWM pwm;
    size_t width = 0;
    double z_score = 0.0;
    double threshold = 0.0;
    double enrichment = -std::numeric_limits<double>::infinity();
    int primary_hits = 0;
    int control_hits = 0;
    int em_iterations = 0;
};

using MaskRow = std::vector<uint8_t>;
using MaskMatrix = std::vector<MaskRow>;

class MarkovOrder3 {
public:
    void train(const SequenceList& sequences, double pseudocount = 1.0) {
        std::array<double, ALPHABET_SIZE> base_counts{};
        base_counts.fill(pseudocount);

        for (auto& row : transition_counts_) {
            row.fill(pseudocount);
        }

        for (const auto& seq : sequences) {
            for (size_t pos = 0; pos < seq.size(); ++pos) {
                const size_t base_idx = base_to_index(seq[pos]);
                if (base_idx == SIZE_MAX) {
                    continue;
                }
                base_counts[base_idx] += 1.0;

                if (pos < 3) {
                    continue;
                }

                const size_t b1 = base_to_index(seq[pos - 3]);
                const size_t b2 = base_to_index(seq[pos - 2]);
                const size_t b3 = base_to_index(seq[pos - 1]);
                if (b1 == SIZE_MAX || b2 == SIZE_MAX || b3 == SIZE_MAX) {
                    continue;
                }

                const size_t context = ((b1 * 4) + b2) * 4 + b3;
                transition_counts_[context][base_idx] += 1.0;
            }
        }

        double base_total = 0.0;
        for (double v : base_counts) {
            base_total += v;
        }
        for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
            base_probs_[b] = base_counts[b] / std::max(base_total, kMinProb);
        }

        for (size_t ctx = 0; ctx < 64; ++ctx) {
            double row_total = 0.0;
            for (double v : transition_counts_[ctx]) {
                row_total += v;
            }
            for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
                transition_probs_[ctx][b] = transition_counts_[ctx][b] / std::max(row_total, kMinProb);
            }
        }
    }

    double conditional_prob(const Sequence& seq, size_t pos) const {
        if (pos >= seq.size()) {
            return 0.25;
        }

        const size_t base_idx = base_to_index(seq[pos]);
        if (base_idx == SIZE_MAX) {
            return 0.25;
        }

        if (pos < 3) {
            return std::max(base_probs_[base_idx], kMinProb);
        }

        const size_t b1 = base_to_index(seq[pos - 3]);
        const size_t b2 = base_to_index(seq[pos - 2]);
        const size_t b3 = base_to_index(seq[pos - 1]);
        if (b1 == SIZE_MAX || b2 == SIZE_MAX || b3 == SIZE_MAX) {
            return std::max(base_probs_[base_idx], kMinProb);
        }

        const size_t context = ((b1 * 4) + b2) * 4 + b3;
        return std::max(transition_probs_[context][base_idx], kMinProb);
    }

private:
    std::array<std::array<double, ALPHABET_SIZE>, 64> transition_counts_{};
    std::array<std::array<double, ALPHABET_SIZE>, 64> transition_probs_{};
    std::array<double, ALPHABET_SIZE> base_probs_ = {0.25, 0.25, 0.25, 0.25};
};

void print_separator(char ch = '=', int width = 90) {
    std::cout << std::string(width, ch) << '\n';
}

double safe_log2(double p) {
    return std::log2(std::max(p, kMinProb));
}

double z_score_two_proportion(int x1, int n1, int x2, int n2) {
    if (n1 <= 0 || n2 <= 0) {
        return 0.0;
    }

    const double p1 = static_cast<double>(x1) / static_cast<double>(n1);
    const double p2 = static_cast<double>(x2) / static_cast<double>(n2);
    const double pooled = static_cast<double>(x1 + x2) / static_cast<double>(n1 + n2);
    const double var = pooled * (1.0 - pooled) * (1.0 / n1 + 1.0 / n2);
    if (var <= 0.0) {
        return 0.0;
    }
    return (p1 - p2) / std::sqrt(var);
}

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

MaskMatrix init_mask(const SequenceList& sequences) {
    MaskMatrix mask;
    mask.reserve(sequences.size());
    for (const auto& seq : sequences) {
        mask.emplace_back(seq.size(), static_cast<uint8_t>(0));
    }
    return mask;
}

SequenceList apply_mask(const SequenceList& sequences, const MaskMatrix& mask) {
    SequenceList masked = sequences;
    const size_t n = std::min(masked.size(), mask.size());
    for (size_t i = 0; i < n; ++i) {
        const size_t m = std::min(masked[i].size(), mask[i].size());
        for (size_t pos = 0; pos < m; ++pos) {
            if (mask[i][pos] != 0) {
                masked[i][pos] = MASK_CHAR;
            }
        }
    }
    return masked;
}

bool window_is_masked(const MaskRow& mask, size_t start, size_t width) {
    if (start + width > mask.size()) {
        return true;
    }
    for (size_t pos = start; pos < start + width; ++pos) {
        if (mask[pos] != 0) {
            return true;
        }
    }
    return false;
}

void mark_window(MaskRow& mask, size_t start, size_t width) {
    if (start + width > mask.size()) {
        return;
    }
    for (size_t pos = start; pos < start + width; ++pos) {
        mask[pos] = 1;
    }
}

double masked_fraction(const MaskMatrix& mask) {
    size_t total = 0;
    size_t masked = 0;
    for (const auto& row : mask) {
        total += row.size();
        for (uint8_t bit : row) {
            masked += (bit != 0) ? 1 : 0;
        }
    }
    return (total == 0) ? 0.0 : static_cast<double>(masked) / static_cast<double>(total);
}

std::vector<std::string> collect_hd_sites(
    const std::string& seed,
    const SequenceList& masked_primary,
    size_t max_hamming = 1
) {
    std::vector<std::string> sites;
    if (seed.empty()) {
        return sites;
    }

    sites.reserve(masked_primary.size());
    for (const auto& seq : masked_primary) {
        if (seq.size() < seed.size()) {
            continue;
        }

        size_t best_hd = std::numeric_limits<size_t>::max();
        size_t best_pos = std::numeric_limits<size_t>::max();

        for (size_t pos = 0; pos + seed.size() <= seq.size(); ++pos) {
            const std::string window = seq.substr(pos, seed.size());
            if (!is_valid_kmer(window)) {
                continue;
            }
            const size_t hd = hamming_distance(window, seed);
            if (hd < best_hd) {
                best_hd = hd;
                best_pos = pos;
                if (hd == 0) {
                    break;
                }
            }
        }

        if (best_pos != std::numeric_limits<size_t>::max() && best_hd <= max_hamming) {
            sites.push_back(seq.substr(best_pos, seed.size()));
        }
    }

    return sites;
}

double score_window_llr(
    const Sequence& seq,
    size_t start,
    const PWM& pwm,
    const MarkovOrder3& background
) {
    if (pwm.width == 0 || start + pwm.width > seq.size()) {
        return -std::numeric_limits<double>::infinity();
    }

    double score = 0.0;
    for (size_t j = 0; j < pwm.width; ++j) {
        const size_t pos = start + j;
        const size_t base_idx = base_to_index(seq[pos]);
        if (base_idx == SIZE_MAX) {
            return -std::numeric_limits<double>::infinity();
        }

        const double pwm_prob = std::max(pwm.matrix[base_idx][j], kMinProb);
        const double bg_prob = std::max(background.conditional_prob(seq, pos), kMinProb);
        score += safe_log2(pwm_prob) - safe_log2(bg_prob);
    }

    return score;
}

MatchResult find_best_site(
    const Sequence& seq,
    const MaskRow& mask,
    const PWM& pwm,
    const MarkovOrder3& background
) {
    MatchResult best;
    if (pwm.width == 0 || seq.size() < pwm.width) {
        return best;
    }

    for (size_t start = 0; start + pwm.width <= seq.size(); ++start) {
        if (window_is_masked(mask, start, pwm.width)) {
            continue;
        }

        const double score = score_window_llr(seq, start, pwm, background);
        if (score > best.score) {
            best.score = score;
            best.position = static_cast<int>(start);
        }
    }

    return best;
}

std::vector<MatchResult> scan_best_sites(
    const SequenceList& sequences,
    const MaskMatrix& mask,
    const PWM& pwm,
    const MarkovOrder3& background
) {
    std::vector<MatchResult> results;
    results.reserve(sequences.size());

    const size_t n = std::min(sequences.size(), mask.size());
    for (size_t i = 0; i < n; ++i) {
        results.push_back(find_best_site(sequences[i], mask[i], pwm, background));
    }
    return results;
}

ThresholdResult choose_best_threshold(
    const std::vector<MatchResult>& primary_matches,
    const std::vector<MatchResult>& control_matches
) {
    std::vector<double> candidates;
    candidates.reserve(primary_matches.size() + control_matches.size());
    for (const auto& m : primary_matches) {
        if (std::isfinite(m.score)) {
            candidates.push_back(m.score);
        }
    }
    for (const auto& m : control_matches) {
        if (std::isfinite(m.score)) {
            candidates.push_back(m.score);
        }
    }

    ThresholdResult best;
    if (candidates.empty()) {
        return best;
    }

    std::sort(candidates.begin(), candidates.end(), std::greater<double>());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    const size_t stride = std::max<size_t>(1, candidates.size() / 128);
    for (size_t i = 0; i < candidates.size(); i += stride) {
        const double threshold = candidates[i];

        int hp = 0;
        int hc = 0;
        for (const auto& m : primary_matches) {
            if (m.score >= threshold) {
                ++hp;
            }
        }
        for (const auto& m : control_matches) {
            if (m.score >= threshold) {
                ++hc;
            }
        }
        if (hp <= 0) {
            continue;
        }

        const double p_rate = (static_cast<double>(hp) + 1.0) /
                              (static_cast<double>(primary_matches.size()) + 2.0);
        const double c_rate = (static_cast<double>(hc) + 1.0) /
                              (static_cast<double>(control_matches.size()) + 2.0);
        const double enrichment = p_rate / std::max(c_rate, kMinProb);

        if (!best.valid || enrichment > best.enrichment ||
            (std::abs(enrichment - best.enrichment) < 1e-9 && hp > best.primary_hits)) {
            best.valid = true;
            best.threshold = threshold;
            best.enrichment = enrichment;
            best.primary_hits = hp;
            best.control_hits = hc;
        }
    }

    return best;
}

std::vector<std::string> collect_sites_for_mstep(
    const SequenceList& sequences,
    const std::vector<MatchResult>& matches,
    double threshold,
    size_t width,
    size_t max_sites
) {
    std::vector<std::string> sites;
    sites.reserve(std::min(sequences.size(), max_sites));

    const size_t n = std::min(sequences.size(), matches.size());
    for (size_t i = 0; i < n; ++i) {
        if (sites.size() >= max_sites) {
            break;
        }

        const auto& m = matches[i];
        if (m.position < 0 || m.score < threshold) {
            continue;
        }

        const size_t pos = static_cast<size_t>(m.position);
        if (pos + width > sequences[i].size()) {
            continue;
        }

        const std::string window = sequences[i].substr(pos, width);
        if (is_valid_kmer(window)) {
            sites.push_back(window);
        }
    }

    return sites;
}

std::vector<SeedCandidate> build_seed_pool(
    const SequenceList& masked_primary,
    const SequenceList& masked_control,
    const PipelineConfig& cfg
) {
    std::vector<SeedCandidate> pool;

    for (size_t k : cfg.k_values) {
        const KmerCounts pos_counts = count_kmers(masked_primary, k);
        const KmerCounts neg_counts = count_kmers(masked_control, k);
        if (pos_counts.empty()) {
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
                static_cast<int>(masked_primary.size()),
                x2,
                static_cast<int>(masked_control.size())
            );

            if (z >= cfg.sub_significant_z) {
                local.push_back({kmer, k, z, x1, x2});
            }
        }

        std::sort(
            local.begin(),
            local.end(),
            [](const SeedCandidate& a, const SeedCandidate& b) {
                if (a.z_score != b.z_score) {
                    return a.z_score > b.z_score;
                }
                return a.kmer < b.kmer;
            }
        );
        if (local.size() > cfg.top_seeds_per_k) {
            local.resize(cfg.top_seeds_per_k);
        }

        pool.insert(pool.end(), local.begin(), local.end());
    }

    std::sort(
        pool.begin(),
        pool.end(),
        [](const SeedCandidate& a, const SeedCandidate& b) {
            if (a.z_score != b.z_score) {
                return a.z_score > b.z_score;
            }
            return a.kmer < b.kmer;
        }
    );

    return pool;
}

std::optional<MotifRoundResult> refine_seed_with_em(
    const SeedCandidate& seed,
    const SequenceList& primary,
    const SequenceList& control,
    const MaskMatrix& primary_mask,
    const MaskMatrix& control_mask,
    const MarkovOrder3& background,
    const PipelineConfig& cfg
) {
    const SequenceList masked_primary = apply_mask(primary, primary_mask);
    const std::vector<std::string> initial_sites = collect_hd_sites(seed.kmer, masked_primary, 1);
    if (initial_sites.size() < cfg.min_sites_for_em) {
        return std::nullopt;
    }

    PWM pwm = build_pwm_from_sites(initial_sites, cfg.pseudocount);
    if (pwm.width == 0) {
        return std::nullopt;
    }

    MotifRoundResult best;
    best.seed = seed.kmer;
    best.z_score = seed.z_score;
    best.width = pwm.width;

    double prev_enrichment = -std::numeric_limits<double>::infinity();
    int stale_steps = 0;

    for (int iter = 0; iter < cfg.max_em_iterations; ++iter) {
        const auto primary_matches = scan_best_sites(primary, primary_mask, pwm, background);
        const auto control_matches = scan_best_sites(control, control_mask, pwm, background);

        const ThresholdResult threshold = choose_best_threshold(primary_matches, control_matches);
        if (!threshold.valid || threshold.primary_hits < static_cast<int>(cfg.min_sites_for_em)) {
            break;
        }

        const auto selected_sites = collect_sites_for_mstep(
            primary,
            primary_matches,
            threshold.threshold,
            pwm.width,
            cfg.max_sites_for_mstep
        );
        if (selected_sites.size() < cfg.min_sites_for_em) {
            break;
        }

        PWM next_pwm = build_pwm_from_sites(selected_sites, cfg.pseudocount);
        if (next_pwm.width == 0) {
            break;
        }

        MotifRoundResult current;
        current.seed = seed.kmer;
        current.z_score = seed.z_score;
        current.pwm = next_pwm;
        current.width = next_pwm.width;
        current.consensus = get_consensus(next_pwm);
        current.threshold = threshold.threshold;
        current.enrichment = threshold.enrichment;
        current.primary_hits = threshold.primary_hits;
        current.control_hits = threshold.control_hits;
        current.em_iterations = iter + 1;

        if (current.enrichment > best.enrichment) {
            best = current;
        }

        if (prev_enrichment > -std::numeric_limits<double>::infinity() &&
            current.enrichment - prev_enrichment < cfg.em_min_improvement) {
            ++stale_steps;
        } else {
            stale_steps = 0;
        }

        prev_enrichment = current.enrichment;
        pwm = std::move(next_pwm);

        if (stale_steps >= cfg.em_patience) {
            break;
        }
    }

    if (best.enrichment < cfg.min_enrichment_to_report) {
        return std::nullopt;
    }
    return best;
}

int erase_by_pwm(
    const SequenceList& sequences,
    MaskMatrix& mask,
    const PWM& pwm,
    const MarkovOrder3& background,
    double threshold
) {
    if (pwm.width == 0) {
        return 0;
    }

    int hit_count = 0;
    const size_t n = std::min(sequences.size(), mask.size());
    for (size_t i = 0; i < n; ++i) {
        const auto& seq = sequences[i];
        if (seq.size() < pwm.width) {
            continue;
        }

        for (size_t start = 0; start + pwm.width <= seq.size(); ++start) {
            if (window_is_masked(mask[i], start, pwm.width)) {
                continue;
            }
            const double score = score_window_llr(seq, start, pwm, background);
            if (score >= threshold) {
                mark_window(mask[i], start, pwm.width);
                ++hit_count;
            }
        }
    }

    return hit_count;
}

std::pair<SequenceList, SequenceList> build_dataset(const PipelineConfig& cfg) {
    if (cfg.use_real_data) {
        return generate_real_test_data(
            cfg.baseline_file,
            cfg.sites_file,
            cfg.num_sequences,
            cfg.segment_length,
            cfg.seed
        );
    }

    DataGenerator gen(cfg.seed);
    return gen.generate_test_data(
        cfg.num_sequences,
        cfg.seq_length,
        cfg.synthetic_motif,
        cfg.synthetic_injection_rate
    );
}

PipelineConfig parse_args(int argc, char* argv[]) {
    PipelineConfig cfg;

    auto parse_size_t = [](const std::string& s) {
        size_t pos = 0;
        const size_t v = std::stoull(s, &pos, 10);
        if (pos != s.size()) {
            throw std::invalid_argument("Invalid integer");
        }
        return v;
    };
    auto parse_int = [](const std::string& s) {
        size_t pos = 0;
        const int v = std::stoi(s, &pos, 10);
        if (pos != s.size()) {
            throw std::invalid_argument("Invalid integer");
        }
        return v;
    };
    auto parse_double = [](const std::string& s) {
        size_t pos = 0;
        const double v = std::stod(s, &pos);
        if (pos != s.size()) {
            throw std::invalid_argument("Invalid float");
        }
        return v;
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                throw std::invalid_argument(std::string("Missing value for ") + name);
            }
            return argv[++i];
        };

        if (arg == "--real") {
            cfg.use_real_data = true;
        } else if (arg == "--baseline") {
            cfg.baseline_file = require_value("--baseline");
        } else if (arg == "--sites") {
            cfg.sites_file = require_value("--sites");
        } else if (arg == "--n") {
            cfg.num_sequences = parse_size_t(require_value("--n"));
        } else if (arg == "--l") {
            cfg.seq_length = parse_size_t(require_value("--l"));
        } else if (arg == "--segment") {
            cfg.segment_length = parse_size_t(require_value("--segment"));
        } else if (arg == "--seed") {
            cfg.seed = static_cast<unsigned int>(parse_int(require_value("--seed")));
        } else if (arg == "--max-motifs") {
            cfg.max_motifs = parse_int(require_value("--max-motifs"));
        } else if (arg == "--motif") {
            cfg.synthetic_motif = require_value("--motif");
        } else if (arg == "--rate") {
            cfg.synthetic_injection_rate = parse_double(require_value("--rate"));
        } else if (arg == "--quick") {
            cfg.num_sequences = 3000;
            cfg.seq_length = 80;
            cfg.max_motifs = 3;
            cfg.max_seed_trials_per_round = 40;
        } else if (arg == "-h" || arg == "--help") {
            std::cout
                << "Usage: motif_discovery [options]\n\n"
                << "Options:\n"
                << "  --real                 Use real DNA dataset mode\n"
                << "  --baseline <file>      Baseline FASTA for real mode\n"
                << "  --sites <file>         Binding sites FASTA for real mode\n"
                << "  --n <int>              Number of sequences per set\n"
                << "  --l <int>              Synthetic sequence length\n"
                << "  --segment <int>        Real-mode chop segment length\n"
                << "  --motif <str>          Synthetic hidden motif\n"
                << "  --rate <float>         Synthetic injection rate [0,1]\n"
                << "  --seed <int>           RNG seed\n"
                << "  --max-motifs <int>     Number of motifs to discover\n"
                << "  --quick                Small quick-run preset\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }

    if (cfg.synthetic_injection_rate < 0.0) {
        cfg.synthetic_injection_rate = 0.0;
    }
    if (cfg.synthetic_injection_rate > 1.0) {
        cfg.synthetic_injection_rate = 1.0;
    }
    cfg.max_motifs = std::max(1, cfg.max_motifs);
    return cfg;
}

} // namespace

int main(int argc, char* argv[]) {
    PipelineConfig cfg;
    try {
        cfg = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << '\n';
        return 2;
    }

    print_separator();
    std::cout << "Integrated Motif Discovery Software Reference\n";
    print_separator();
    std::cout << "Mode                : " << (cfg.use_real_data ? "real" : "synthetic") << '\n';
    std::cout << "Sequences per set   : " << cfg.num_sequences << '\n';
    if (!cfg.use_real_data) {
        std::cout << "Sequence length      : " << cfg.seq_length << '\n';
        std::cout << "Injected motif       : " << cfg.synthetic_motif << '\n';
    } else {
        std::cout << "Segment length       : " << cfg.segment_length << '\n';
        std::cout << "Baseline file        : " << cfg.baseline_file << '\n';
        std::cout << "Sites file           : " << cfg.sites_file << '\n';
    }
    std::cout << "k values             : ";
    for (size_t i = 0; i < cfg.k_values.size(); ++i) {
        std::cout << cfg.k_values[i] << (i + 1 == cfg.k_values.size() ? '\n' : ',');
    }
    std::cout << "z thresholds         : sig>=" << cfg.significant_z
              << ", sub>=" << cfg.sub_significant_z << '\n';
    std::cout << "max motifs           : " << cfg.max_motifs << "\n\n";

    SequenceList primary;
    SequenceList control;
    try {
        std::tie(primary, control) = build_dataset(cfg);
    } catch (const std::exception& e) {
        std::cerr << "Dataset error: " << e.what() << '\n';
        return 1;
    }

    std::cout << "Loaded dataset: primary=" << primary.size() << ", control=" << control.size() << "\n";
    if (primary.empty() || control.empty()) {
        std::cerr << "Empty dataset; aborting.\n";
        return 1;
    }

    MarkovOrder3 background;
    background.train(control);

    MaskMatrix primary_mask = init_mask(primary);
    MaskMatrix control_mask = init_mask(control);

    std::vector<MotifRoundResult> discovered;
    discovered.reserve(static_cast<size_t>(cfg.max_motifs));

    for (int round = 0; round < cfg.max_motifs; ++round) {
        print_separator('-');
        std::cout << "Round " << (round + 1) << '\n';
        print_separator('-');

        const SequenceList masked_primary = apply_mask(primary, primary_mask);
        const SequenceList masked_control = apply_mask(control, control_mask);

        std::vector<SeedCandidate> seed_pool = build_seed_pool(masked_primary, masked_control, cfg);
        if (seed_pool.empty()) {
            std::cout << "No candidate seed passed z-thresholds.\n";
            break;
        }

        std::optional<MotifRoundResult> best;
        size_t tried = 0;
        for (const auto& seed : seed_pool) {
            if (tried++ >= cfg.max_seed_trials_per_round) {
                break;
            }

            const auto refined = refine_seed_with_em(
                seed,
                primary,
                control,
                primary_mask,
                control_mask,
                background,
                cfg
            );

            if (!refined.has_value()) {
                continue;
            }

            if (!best.has_value() || refined->enrichment > best->enrichment) {
                best = refined;
            }
        }

        if (!best.has_value()) {
            std::cout << "No seed converged to a valid motif in this round.\n";
            break;
        }

        const int erased_primary = erase_by_pwm(primary, primary_mask, best->pwm, background, best->threshold);
        const int erased_control = erase_by_pwm(control, control_mask, best->pwm, background, best->threshold);

        std::cout << "Seed                : " << best->seed << " (z=" << std::fixed << std::setprecision(3)
                  << best->z_score << ")\n";
        std::cout << "Consensus           : " << best->consensus << " (W=" << best->width << ")\n";
        std::cout << "Enrichment          : " << std::fixed << std::setprecision(3) << best->enrichment << '\n';
        std::cout << "Threshold (LLR)     : " << std::fixed << std::setprecision(3) << best->threshold << '\n';
        std::cout << "Hits (P/C)          : " << best->primary_hits << " / " << best->control_hits << '\n';
        std::cout << "Erased windows (P/C): " << erased_primary << " / " << erased_control << '\n';
        std::cout << "EM iterations       : " << best->em_iterations << '\n';

        discovered.push_back(*best);

        if (erased_primary == 0) {
            std::cout << "No new primary windows erased; stopping.\n";
            break;
        }

        const double frac = masked_fraction(primary_mask);
        std::cout << "Primary masked frac : " << std::fixed << std::setprecision(2) << (frac * 100.0) << "%\n";
        if (frac >= cfg.max_mask_fraction) {
            std::cout << "Mask fraction limit reached; stopping.\n";
            break;
        }
    }

    print_separator();
    std::cout << "Discovered motifs: " << discovered.size() << '\n';
    print_separator();
    for (size_t i = 0; i < discovered.size(); ++i) {
        const auto& m = discovered[i];
        std::cout << (i + 1) << '\t' << m.consensus << '\t'
                  << "W=" << m.width << '\t'
                  << "z=" << std::fixed << std::setprecision(3) << m.z_score << '\t'
                  << "enrich=" << std::fixed << std::setprecision(3) << m.enrichment << '\n';
    }

    return 0;
}
