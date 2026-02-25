/**
 * @file benchmark.cpp
 * @brief Integrated pipeline benchmark for software validation before FPGA design.
 */

#include "motif.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>

using namespace motif;
using Clock = std::chrono::high_resolution_clock;

namespace {

constexpr double kMinProb = 1e-12;

struct BenchmarkConfig {
    bool real_mode = false;
    std::string baseline_file = "data/upstream5000.fa";
    std::string sites_file = "data/MA0007.2.sites";
    size_t segment_length = 1000;

    std::vector<std::pair<size_t, size_t>> synthetic_shapes = {
        {3000, 80},
        {6000, 100},
        {12000, 100},
        {24000, 120}
    };

    int runs = 3;
    unsigned int seed = 42;
    std::string motif = "TATATA";
    double rate = 0.45;
};

struct SeedCandidate {
    std::string kmer;
    size_t width = 0;
    double z_score = 0.0;
};

struct MatchResult {
    int pos = -1;
    double score = -std::numeric_limits<double>::infinity();
};

struct ThresholdResult {
    bool valid = false;
    double threshold = 0.0;
    double enrichment = 0.0;
    int primary_hits = 0;
    int control_hits = 0;
};

struct MotifSummary {
    std::string consensus;
    double enrichment = 0.0;
};

struct PipelineConfig {
    std::vector<size_t> k_values = {6, 8, 9, 10, 12};
    double sub_significant_z = 3.5;
    size_t top_seeds_per_k = 12;
    size_t max_seed_trials = 24;
    int max_motifs = 3;
    int max_em_iterations = 10;
    int em_patience = 2;
    double em_min_improvement = 1e-3;
    size_t min_sites_for_em = 6;
    size_t max_sites_for_mstep = 20000;
    double pseudocount = 0.2;
    double min_enrichment = 1.2;
    double max_mask_fraction = 0.85;
};

using MaskRow = std::vector<uint8_t>;
using MaskMatrix = std::vector<MaskRow>;

class MarkovOrder3 {
public:
    void train(const SequenceList& sequences) {
        std::array<double, ALPHABET_SIZE> base_counts = {1.0, 1.0, 1.0, 1.0};
        for (auto& row : trans_counts_) {
            row = {1.0, 1.0, 1.0, 1.0};
        }

        for (const auto& seq : sequences) {
            for (size_t i = 0; i < seq.size(); ++i) {
                const size_t b = base_to_index(seq[i]);
                if (b == SIZE_MAX) {
                    continue;
                }
                base_counts[b] += 1.0;

                if (i < 3) {
                    continue;
                }

                const size_t b1 = base_to_index(seq[i - 3]);
                const size_t b2 = base_to_index(seq[i - 2]);
                const size_t b3 = base_to_index(seq[i - 1]);
                if (b1 == SIZE_MAX || b2 == SIZE_MAX || b3 == SIZE_MAX) {
                    continue;
                }
                const size_t ctx = ((b1 * 4) + b2) * 4 + b3;
                trans_counts_[ctx][b] += 1.0;
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
            double sum = 0.0;
            for (double v : trans_counts_[ctx]) {
                sum += v;
            }
            for (size_t b = 0; b < ALPHABET_SIZE; ++b) {
                trans_probs_[ctx][b] = trans_counts_[ctx][b] / std::max(sum, kMinProb);
            }
        }
    }

    double conditional_prob(const Sequence& seq, size_t pos) const {
        const size_t b = (pos < seq.size()) ? base_to_index(seq[pos]) : SIZE_MAX;
        if (b == SIZE_MAX) {
            return 0.25;
        }
        if (pos < 3) {
            return std::max(base_probs_[b], kMinProb);
        }

        const size_t b1 = base_to_index(seq[pos - 3]);
        const size_t b2 = base_to_index(seq[pos - 2]);
        const size_t b3 = base_to_index(seq[pos - 1]);
        if (b1 == SIZE_MAX || b2 == SIZE_MAX || b3 == SIZE_MAX) {
            return std::max(base_probs_[b], kMinProb);
        }
        const size_t ctx = ((b1 * 4) + b2) * 4 + b3;
        return std::max(trans_probs_[ctx][b], kMinProb);
    }

private:
    std::array<std::array<double, ALPHABET_SIZE>, 64> trans_counts_{};
    std::array<std::array<double, ALPHABET_SIZE>, 64> trans_probs_{};
    std::array<double, ALPHABET_SIZE> base_probs_ = {0.25, 0.25, 0.25, 0.25};
};

double z_score_two_proportion(int x1, int n1, int x2, int n2) {
    if (n1 <= 0 || n2 <= 0) {
        return 0.0;
    }
    const double p1 = static_cast<double>(x1) / n1;
    const double p2 = static_cast<double>(x2) / n2;
    const double pooled = static_cast<double>(x1 + x2) / (n1 + n2);
    const double var = pooled * (1.0 - pooled) * (1.0 / n1 + 1.0 / n2);
    return (var <= 0.0) ? 0.0 : (p1 - p2) / std::sqrt(var);
}

MaskMatrix init_mask(const SequenceList& seqs) {
    MaskMatrix mask;
    mask.reserve(seqs.size());
    for (const auto& seq : seqs) {
        mask.emplace_back(seq.size(), static_cast<uint8_t>(0));
    }
    return mask;
}

SequenceList apply_mask(const SequenceList& seqs, const MaskMatrix& mask) {
    SequenceList out = seqs;
    for (size_t i = 0; i < std::min(out.size(), mask.size()); ++i) {
        for (size_t j = 0; j < std::min(out[i].size(), mask[i].size()); ++j) {
            if (mask[i][j] != 0) {
                out[i][j] = MASK_CHAR;
            }
        }
    }
    return out;
}

bool masked_window(const MaskRow& mask, size_t start, size_t width) {
    if (start + width > mask.size()) {
        return true;
    }
    for (size_t i = start; i < start + width; ++i) {
        if (mask[i] != 0) {
            return true;
        }
    }
    return false;
}

void mark_window(MaskRow& mask, size_t start, size_t width) {
    if (start + width > mask.size()) {
        return;
    }
    for (size_t i = start; i < start + width; ++i) {
        mask[i] = 1;
    }
}

double score_window_llr(const Sequence& seq, size_t start, const PWM& pwm, const MarkovOrder3& bg) {
    if (pwm.width == 0 || start + pwm.width > seq.size()) {
        return -std::numeric_limits<double>::infinity();
    }

    double score = 0.0;
    for (size_t j = 0; j < pwm.width; ++j) {
        const size_t base_idx = base_to_index(seq[start + j]);
        if (base_idx == SIZE_MAX) {
            return -std::numeric_limits<double>::infinity();
        }
        const double pwm_p = std::max(pwm.matrix[base_idx][j], kMinProb);
        const double bg_p = std::max(bg.conditional_prob(seq, start + j), kMinProb);
        score += std::log2(pwm_p) - std::log2(bg_p);
    }
    return score;
}

MatchResult best_site(const Sequence& seq, const MaskRow& mask, const PWM& pwm, const MarkovOrder3& bg) {
    MatchResult best;
    if (seq.size() < pwm.width || pwm.width == 0) {
        return best;
    }

    for (size_t start = 0; start + pwm.width <= seq.size(); ++start) {
        if (masked_window(mask, start, pwm.width)) {
            continue;
        }
        const double s = score_window_llr(seq, start, pwm, bg);
        if (s > best.score) {
            best.score = s;
            best.pos = static_cast<int>(start);
        }
    }
    return best;
}

std::vector<MatchResult> scan_set(const SequenceList& seqs, const MaskMatrix& mask, const PWM& pwm, const MarkovOrder3& bg) {
    std::vector<MatchResult> out;
    out.reserve(seqs.size());
    for (size_t i = 0; i < std::min(seqs.size(), mask.size()); ++i) {
        out.push_back(best_site(seqs[i], mask[i], pwm, bg));
    }
    return out;
}

ThresholdResult choose_threshold(const std::vector<MatchResult>& pm, const std::vector<MatchResult>& cm) {
    std::vector<double> cand;
    cand.reserve(pm.size() + cm.size());
    for (const auto& m : pm) {
        if (std::isfinite(m.score)) cand.push_back(m.score);
    }
    for (const auto& m : cm) {
        if (std::isfinite(m.score)) cand.push_back(m.score);
    }
    std::sort(cand.begin(), cand.end(), std::greater<double>());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

    ThresholdResult best;
    if (cand.empty()) {
        return best;
    }

    const size_t stride = std::max<size_t>(1, cand.size() / 128);
    for (size_t i = 0; i < cand.size(); i += stride) {
        const double t = cand[i];
        int hp = 0;
        int hc = 0;
        for (const auto& m : pm) if (m.score >= t) ++hp;
        for (const auto& m : cm) if (m.score >= t) ++hc;
        if (hp == 0) continue;

        const double p_rate = (hp + 1.0) / (pm.size() + 2.0);
        const double c_rate = (hc + 1.0) / (cm.size() + 2.0);
        const double enrich = p_rate / std::max(c_rate, kMinProb);

        if (!best.valid || enrich > best.enrichment ||
            (std::abs(enrich - best.enrichment) < 1e-9 && hp > best.primary_hits)) {
            best.valid = true;
            best.threshold = t;
            best.enrichment = enrich;
            best.primary_hits = hp;
            best.control_hits = hc;
        }
    }
    return best;
}

size_t hamming_distance(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return std::numeric_limits<size_t>::max();
    size_t d = 0;
    for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) ++d;
    return d;
}

std::vector<std::string> collect_hd_sites(const std::string& seed, const SequenceList& masked_primary) {
    std::vector<std::string> sites;
    sites.reserve(masked_primary.size());
    for (const auto& seq : masked_primary) {
        if (seq.size() < seed.size()) continue;
        size_t best_hd = std::numeric_limits<size_t>::max();
        size_t best_pos = std::numeric_limits<size_t>::max();
        for (size_t p = 0; p + seed.size() <= seq.size(); ++p) {
            const std::string w = seq.substr(p, seed.size());
            if (!is_valid_kmer(w)) continue;
            const size_t hd = hamming_distance(w, seed);
            if (hd < best_hd) {
                best_hd = hd;
                best_pos = p;
                if (hd == 0) break;
            }
        }
        if (best_pos != std::numeric_limits<size_t>::max() && best_hd <= 1) {
            sites.push_back(seq.substr(best_pos, seed.size()));
        }
    }
    return sites;
}

std::vector<SeedCandidate> build_seed_pool(const SequenceList& masked_primary, const SequenceList& masked_control, const PipelineConfig& cfg) {
    std::vector<SeedCandidate> pool;

    for (size_t k : cfg.k_values) {
        const auto pc = count_kmers(masked_primary, k);
        const auto nc = count_kmers(masked_control, k);
        std::set<std::string> all;
        for (const auto& [kmer, _] : pc) all.insert(kmer);
        for (const auto& [kmer, _] : nc) all.insert(kmer);

        std::vector<SeedCandidate> local;
        for (const auto& kmer : all) {
            const int x1 = pc.count(kmer) ? pc.at(kmer) : 0;
            const int x2 = nc.count(kmer) ? nc.at(kmer) : 0;
            const double z = z_score_two_proportion(x1, static_cast<int>(masked_primary.size()), x2, static_cast<int>(masked_control.size()));
            if (z >= cfg.sub_significant_z) {
                local.push_back({kmer, k, z});
            }
        }
        std::sort(local.begin(), local.end(), [](const SeedCandidate& a, const SeedCandidate& b) {
            if (a.z_score != b.z_score) return a.z_score > b.z_score;
            return a.kmer < b.kmer;
        });
        if (local.size() > cfg.top_seeds_per_k) local.resize(cfg.top_seeds_per_k);
        pool.insert(pool.end(), local.begin(), local.end());
    }

    std::sort(pool.begin(), pool.end(), [](const SeedCandidate& a, const SeedCandidate& b) {
        if (a.z_score != b.z_score) return a.z_score > b.z_score;
        return a.kmer < b.kmer;
    });
    return pool;
}

std::optional<MotifSummary> refine_seed(
    const SeedCandidate& seed,
    const SequenceList& primary,
    const SequenceList& control,
    const MaskMatrix& pmask,
    const MaskMatrix& cmask,
    const MarkovOrder3& bg,
    const PipelineConfig& cfg,
    PWM* out_pwm,
    double* out_threshold
) {
    const SequenceList masked_primary = apply_mask(primary, pmask);
    const auto init_sites = collect_hd_sites(seed.kmer, masked_primary);
    if (init_sites.size() < cfg.min_sites_for_em) return std::nullopt;

    PWM pwm = build_pwm_from_sites(init_sites, cfg.pseudocount);
    if (pwm.width == 0) return std::nullopt;

    MotifSummary best;
    bool found = false;
    double best_thr = 0.0;
    double prev_enr = -std::numeric_limits<double>::infinity();
    int stale = 0;

    for (int iter = 0; iter < cfg.max_em_iterations; ++iter) {
        const auto pbest = scan_set(primary, pmask, pwm, bg);
        const auto cbest = scan_set(control, cmask, pwm, bg);
        const auto thr = choose_threshold(pbest, cbest);
        if (!thr.valid || thr.primary_hits < static_cast<int>(cfg.min_sites_for_em)) break;

        std::vector<std::string> sites;
        sites.reserve(std::min(primary.size(), cfg.max_sites_for_mstep));
        for (size_t i = 0; i < std::min(primary.size(), pbest.size()) && sites.size() < cfg.max_sites_for_mstep; ++i) {
            if (pbest[i].pos < 0 || pbest[i].score < thr.threshold) continue;
            const size_t pos = static_cast<size_t>(pbest[i].pos);
            if (pos + pwm.width > primary[i].size()) continue;
            const auto w = primary[i].substr(pos, pwm.width);
            if (is_valid_kmer(w)) sites.push_back(w);
        }
        if (sites.size() < cfg.min_sites_for_em) break;

        PWM next = build_pwm_from_sites(sites, cfg.pseudocount);
        if (next.width == 0) break;

        if (!found || thr.enrichment > best.enrichment) {
            best.consensus = get_consensus(next);
            best.enrichment = thr.enrichment;
            best_thr = thr.threshold;
            *out_pwm = next;
            found = true;
        }

        if (prev_enr > -std::numeric_limits<double>::infinity() && thr.enrichment - prev_enr < cfg.em_min_improvement) {
            ++stale;
        } else {
            stale = 0;
        }
        prev_enr = thr.enrichment;
        pwm = std::move(next);
        if (stale >= cfg.em_patience) break;
    }

    if (!found || best.enrichment < cfg.min_enrichment) return std::nullopt;
    *out_threshold = best_thr;
    return best;
}

int erase_hits(const SequenceList& seqs, MaskMatrix& mask, const PWM& pwm, const MarkovOrder3& bg, double thr) {
    int hits = 0;
    for (size_t i = 0; i < std::min(seqs.size(), mask.size()); ++i) {
        const auto& seq = seqs[i];
        if (seq.size() < pwm.width) continue;
        for (size_t s = 0; s + pwm.width <= seq.size(); ++s) {
            if (masked_window(mask[i], s, pwm.width)) continue;
            if (score_window_llr(seq, s, pwm, bg) >= thr) {
                mark_window(mask[i], s, pwm.width);
                ++hits;
            }
        }
    }
    return hits;
}

double masked_fraction(const MaskMatrix& mask) {
    size_t total = 0;
    size_t masked = 0;
    for (const auto& row : mask) {
        total += row.size();
        for (uint8_t b : row) masked += (b != 0);
    }
    return (total == 0) ? 0.0 : static_cast<double>(masked) / total;
}

std::vector<MotifSummary> run_integrated_pipeline(
    const SequenceList& primary,
    const SequenceList& control,
    const PipelineConfig& cfg
) {
    std::vector<MotifSummary> motifs;
    if (primary.empty() || control.empty()) {
        return motifs;
    }

    MarkovOrder3 bg;
    bg.train(control);

    MaskMatrix pmask = init_mask(primary);
    MaskMatrix cmask = init_mask(control);

    for (int round = 0; round < cfg.max_motifs; ++round) {
        const auto mp = apply_mask(primary, pmask);
        const auto mc = apply_mask(control, cmask);
        auto seeds = build_seed_pool(mp, mc, cfg);
        if (seeds.empty()) break;

        std::optional<MotifSummary> best;
        PWM best_pwm;
        double best_thr = 0.0;

        size_t tried = 0;
        for (const auto& seed : seeds) {
            if (tried++ >= cfg.max_seed_trials) break;

            PWM cand_pwm;
            double cand_thr = 0.0;
            auto candidate = refine_seed(seed, primary, control, pmask, cmask, bg, cfg, &cand_pwm, &cand_thr);
            if (!candidate.has_value()) continue;

            if (!best.has_value() || candidate->enrichment > best->enrichment) {
                best = candidate;
                best_pwm = cand_pwm;
                best_thr = cand_thr;
            }
        }

        if (!best.has_value()) break;

        const int erased = erase_hits(primary, pmask, best_pwm, bg, best_thr);
        (void)erase_hits(control, cmask, best_pwm, bg, best_thr);
        motifs.push_back(*best);

        if (erased == 0 || masked_fraction(pmask) >= cfg.max_mask_fraction) {
            break;
        }
    }

    return motifs;
}

struct RunStats {
    double mean_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
};

RunStats summarize(const std::vector<double>& samples) {
    RunStats s;
    if (samples.empty()) return s;
    s.min_ms = *std::min_element(samples.begin(), samples.end());
    s.max_ms = *std::max_element(samples.begin(), samples.end());
    double sum = 0.0;
    for (double v : samples) sum += v;
    s.mean_ms = sum / samples.size();
    return s;
}

void print_separator(char ch = '=', int width = 92) {
    std::cout << std::string(width, ch) << '\n';
}

BenchmarkConfig parse_args(int argc, char* argv[]) {
    BenchmarkConfig cfg;

    auto parse_size_t = [](const std::string& s) {
        size_t p = 0;
        size_t v = std::stoull(s, &p, 10);
        if (p != s.size()) throw std::invalid_argument("invalid integer");
        return v;
    };
    auto parse_int = [](const std::string& s) {
        size_t p = 0;
        int v = std::stoi(s, &p, 10);
        if (p != s.size()) throw std::invalid_argument("invalid integer");
        return v;
    };
    auto parse_double = [](const std::string& s) {
        size_t p = 0;
        double v = std::stod(s, &p);
        if (p != s.size()) throw std::invalid_argument("invalid float");
        return v;
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[++i];
        };

        if (arg == "--quick" || arg == "-q") {
            cfg.synthetic_shapes = {{1000, 80}};
            cfg.runs = 1;
        } else if (arg == "--large") {
            cfg.synthetic_shapes = {{12000, 100}, {24000, 120}, {50000, 120}};
            cfg.runs = 1;
        } else if (arg == "--real") {
            cfg.real_mode = true;
        } else if (arg == "--baseline") {
            cfg.baseline_file = require_value("--baseline");
        } else if (arg == "--sites") {
            cfg.sites_file = require_value("--sites");
        } else if (arg == "--segment") {
            cfg.segment_length = parse_size_t(require_value("--segment"));
        } else if (arg == "--n") {
            const size_t n = parse_size_t(require_value("--n"));
            cfg.synthetic_shapes = {{n, 100}};
        } else if (arg == "--l") {
            const size_t l = parse_size_t(require_value("--l"));
            if (cfg.synthetic_shapes.empty()) cfg.synthetic_shapes.push_back({5000, l});
            for (auto& shape : cfg.synthetic_shapes) shape.second = l;
        } else if (arg == "--runs") {
            cfg.runs = std::max(1, parse_int(require_value("--runs")));
        } else if (arg == "--seed") {
            cfg.seed = static_cast<unsigned int>(std::max(0, parse_int(require_value("--seed"))));
        } else if (arg == "--motif") {
            cfg.motif = require_value("--motif");
        } else if (arg == "--rate") {
            cfg.rate = parse_double(require_value("--rate"));
            if (cfg.rate < 0.0) cfg.rate = 0.0;
            if (cfg.rate > 1.0) cfg.rate = 1.0;
        } else if (arg == "-h" || arg == "--help") {
            std::cout
                << "Usage: motif_benchmark [options]\n\n"
                << "Options:\n"
                << "  -q, --quick          Quick synthetic benchmark\n"
                << "  --large              Large synthetic benchmark\n"
                << "  --real               Real data benchmark\n"
                << "  --baseline <file>    Baseline FASTA file\n"
                << "  --sites <file>       Binding-sites FASTA file\n"
                << "  --segment <int>      Segment length for real data\n"
                << "  --n <int>            Synthetic sequence count\n"
                << "  --l <int>            Synthetic sequence length\n"
                << "  --runs <int>         Number of repeated runs\n"
                << "  --seed <int>         RNG seed\n"
                << "  --motif <str>        Synthetic hidden motif\n"
                << "  --rate <float>       Injection rate [0,1]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }

    return cfg;
}

} // namespace

int main(int argc, char* argv[]) {
    BenchmarkConfig cfg;
    try {
        cfg = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << '\n';
        return 2;
    }

    PipelineConfig pipeline_cfg;

    print_separator();
    std::cout << "Integrated Pipeline Benchmark\n";
    print_separator();
    std::cout << "Mode: " << (cfg.real_mode ? "real" : "synthetic")
              << ", runs per case: " << cfg.runs << "\n\n";

    std::cout << std::setw(12) << "Sequences"
              << std::setw(10) << "Len"
              << std::setw(14) << "Mean(ms)"
              << std::setw(14) << "Min(ms)"
              << std::setw(14) << "Max(ms)"
              << std::setw(10) << "Motifs"
              << std::setw(18) << "TopConsensus"
              << std::setw(12) << "TopEnrich"
              << '\n';
    print_separator('-');

    if (cfg.real_mode) {
        const std::vector<size_t> Ns = {4096, 8192, 16384};
        for (size_t n : Ns) {
            SequenceList primary;
            SequenceList control;
            try {
                std::tie(primary, control) = generate_real_test_data(
                    cfg.baseline_file,
                    cfg.sites_file,
                    n,
                    cfg.segment_length,
                    cfg.seed
                );
            } catch (const std::exception& e) {
                std::cerr << "Real data load error: " << e.what() << '\n';
                return 1;
            }

            std::vector<double> times;
            std::vector<MotifSummary> best_run;
            for (int run = 0; run < cfg.runs; ++run) {
                const auto t0 = Clock::now();
                auto motifs = run_integrated_pipeline(primary, control, pipeline_cfg);
                const auto t1 = Clock::now();
                const double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
                times.push_back(ms);
                if (best_run.empty() || motifs.size() > best_run.size()) {
                    best_run = std::move(motifs);
                }
            }

            const RunStats s = summarize(times);
            const std::string top_consensus = best_run.empty() ? "-" : best_run.front().consensus;
            const double top_enrich = best_run.empty() ? 0.0 : best_run.front().enrichment;

            std::cout << std::setw(12) << n
                      << std::setw(10) << cfg.segment_length
                      << std::setw(14) << std::fixed << std::setprecision(1) << s.mean_ms
                      << std::setw(14) << std::fixed << std::setprecision(1) << s.min_ms
                      << std::setw(14) << std::fixed << std::setprecision(1) << s.max_ms
                      << std::setw(10) << best_run.size()
                      << std::setw(18) << top_consensus
                      << std::setw(12) << std::fixed << std::setprecision(2) << top_enrich
                      << '\n';
        }
    } else {
        DataGenerator gen(cfg.seed);
        for (const auto& [num_seq, seq_len] : cfg.synthetic_shapes) {
            auto [primary, control] = gen.generate_test_data(num_seq, seq_len, cfg.motif, cfg.rate);

            std::vector<double> times;
            std::vector<MotifSummary> best_run;
            for (int run = 0; run < cfg.runs; ++run) {
                const auto t0 = Clock::now();
                auto motifs = run_integrated_pipeline(primary, control, pipeline_cfg);
                const auto t1 = Clock::now();
                const double ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
                times.push_back(ms);
                if (best_run.empty() || motifs.size() > best_run.size()) {
                    best_run = std::move(motifs);
                }
            }

            const RunStats s = summarize(times);
            const std::string top_consensus = best_run.empty() ? "-" : best_run.front().consensus;
            const double top_enrich = best_run.empty() ? 0.0 : best_run.front().enrichment;

            std::cout << std::setw(12) << num_seq
                      << std::setw(10) << seq_len
                      << std::setw(14) << std::fixed << std::setprecision(1) << s.mean_ms
                      << std::setw(14) << std::fixed << std::setprecision(1) << s.min_ms
                      << std::setw(14) << std::fixed << std::setprecision(1) << s.max_ms
                      << std::setw(10) << best_run.size()
                      << std::setw(18) << top_consensus
                      << std::setw(12) << std::fixed << std::setprecision(2) << top_enrich
                      << '\n';
        }
    }

    print_separator();
    std::cout << "Benchmark complete.\n";
    return 0;
}
