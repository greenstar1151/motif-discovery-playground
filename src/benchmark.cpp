/**
 * @file benchmark.cpp
 * @brief Integrated pipeline benchmark — thin wrapper around motif_core.
 */

#include "motif.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

using namespace motif;
using Clock = std::chrono::high_resolution_clock;

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// Configuration (benchmark-specific)
// ─────────────────────────────────────────────────────────────────────────────

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

struct LocalPipelineConfig {
    std::vector<size_t> k_values = {6, 8, 9, 10, 12};
    double sub_significant_z = 3.5;
    size_t top_seeds_per_k = 12;
    size_t max_seed_trials = 24;
    int max_motifs = 3;
    double min_enrichment = 1.2;
    double max_mask_fraction = 0.85;
    EMConfig em{10, 2, 1e-3, 6, 20000, 0.2};
};

struct MotifSummary {
    std::string consensus;
    double enrichment = 0.0;
};

struct RunStats {
    double mean_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

RunStats summarize(const std::vector<double>& samples) {
    RunStats s;
    if (samples.empty()) return s;
    s.min_ms = *std::min_element(samples.begin(), samples.end());
    s.max_ms = *std::max_element(samples.begin(), samples.end());
    double sum = 0.0;
    for (double v : samples) sum += v;
    s.mean_ms = sum / static_cast<double>(samples.size());
    return s;
}

void print_separator(char ch = '=', int width = 92) {
    std::cout << std::string(width, ch) << '\n';
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-seed refinement (delegates to library EM)
// ─────────────────────────────────────────────────────────────────────────────

struct SeedRefineResult {
    MotifSummary summary;
    PWM pwm;
    double threshold = 0.0;
};

std::optional<SeedRefineResult> refine_seed(
    const SeedCandidate& seed,
    const SequenceList& primary,
    const SequenceList& control,
    const MaskMatrix& pmask,
    const MaskMatrix& cmask,
    const MarkovOrder3& bg,
    const LocalPipelineConfig& cfg
) {
    const SequenceList masked_primary = apply_mask(primary, pmask);
    const auto init_sites = collect_hd_sites(seed.kmer, masked_primary);
    if (init_sites.size() < cfg.em.min_sites) {
        return std::nullopt;
    }

    PWM init_pwm = build_pwm_from_sites(init_sites, cfg.em.pseudocount);
    if (init_pwm.width == 0) {
        return std::nullopt;
    }

    auto refined = run_em_refinement(
        init_pwm, primary, control, pmask, cmask, bg, cfg.em, seed.kmer
    );

    if (!refined.has_value() || refined->enrichment < cfg.min_enrichment) {
        return std::nullopt;
    }

    SeedRefineResult result;
    result.summary.consensus = get_consensus(refined->pwm);
    result.summary.enrichment = refined->enrichment;
    result.pwm = refined->pwm;
    result.threshold = refined->learned_threshold;
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Pipeline (multi-round with erasing)
// ─────────────────────────────────────────────────────────────────────────────

std::vector<MotifSummary> run_integrated_pipeline(
    const SequenceList& primary,
    const SequenceList& control,
    const LocalPipelineConfig& cfg
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

        auto seeds = build_seed_pool(
            mp, mc, cfg.k_values, cfg.sub_significant_z, cfg.top_seeds_per_k
        );
        if (seeds.empty()) break;

        std::optional<SeedRefineResult> best;
        size_t tried = 0;
        for (const auto& seed : seeds) {
            if (tried++ >= cfg.max_seed_trials) break;

            auto candidate = refine_seed(
                seed, primary, control, pmask, cmask, bg, cfg
            );
            if (!candidate.has_value()) continue;

            if (!best.has_value() ||
                candidate->summary.enrichment > best->summary.enrichment) {
                best = candidate;
            }
        }

        if (!best.has_value()) break;

        const int erased = erase_by_pwm(
            primary, pmask, best->pwm, bg, best->threshold
        );
        (void)erase_by_pwm(
            control, cmask, best->pwm, bg, best->threshold
        );
        motifs.push_back(best->summary);

        if (erased == 0 || masked_fraction(pmask) >= cfg.max_mask_fraction) {
            break;
        }
    }

    return motifs;
}

// ─────────────────────────────────────────────────────────────────────────────
// Argument parsing
// ─────────────────────────────────────────────────────────────────────────────

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
            if (i + 1 >= argc)
                throw std::invalid_argument(std::string("missing value for ") + name);
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
            if (cfg.synthetic_shapes.empty())
                cfg.synthetic_shapes.push_back({5000, l});
            for (auto& shape : cfg.synthetic_shapes) shape.second = l;
        } else if (arg == "--runs") {
            cfg.runs = std::max(1, parse_int(require_value("--runs")));
        } else if (arg == "--seed") {
            cfg.seed = static_cast<unsigned int>(
                std::max(0, parse_int(require_value("--seed")))
            );
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

// =============================================================================
// main
// =============================================================================

int main(int argc, char* argv[]) {
    BenchmarkConfig cfg;
    try {
        cfg = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << '\n';
        return 2;
    }

    LocalPipelineConfig pipeline_cfg;

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
                    cfg.baseline_file, cfg.sites_file,
                    n, cfg.segment_length, cfg.seed
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
                const double ms =
                    std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
                times.push_back(ms);
                if (best_run.empty() || motifs.size() > best_run.size()) {
                    best_run = std::move(motifs);
                }
            }

            const RunStats s = summarize(times);
            const std::string top_consensus =
                best_run.empty() ? "-" : best_run.front().consensus;
            const double top_enrich =
                best_run.empty() ? 0.0 : best_run.front().enrichment;

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
            auto [primary, control] = gen.generate_test_data(
                num_seq, seq_len, cfg.motif, cfg.rate
            );

            std::vector<double> times;
            std::vector<MotifSummary> best_run;
            for (int run = 0; run < cfg.runs; ++run) {
                const auto t0 = Clock::now();
                auto motifs = run_integrated_pipeline(primary, control, pipeline_cfg);
                const auto t1 = Clock::now();
                const double ms =
                    std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
                times.push_back(ms);
                if (best_run.empty() || motifs.size() > best_run.size()) {
                    best_run = std::move(motifs);
                }
            }

            const RunStats s = summarize(times);
            const std::string top_consensus =
                best_run.empty() ? "-" : best_run.front().consensus;
            const double top_enrich =
                best_run.empty() ? 0.0 : best_run.front().enrichment;

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
