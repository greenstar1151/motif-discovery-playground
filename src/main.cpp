/**
 * @file main.cpp
 * @brief Integrated software reference pipeline for motif discovery.
 *
 * This is a thin CLI wrapper around the shared motif_core library.
 * Pipeline: multi-k z-test seeds → HD≤1 init → EM refinement → erasing loop.
 */

#include "motif.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>

using namespace motif;

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

struct PipelineConfig {
    bool use_real_data = false;
    std::string baseline_file = "data/upstream5000.fa";
    std::string sites_file   = "data/MA0007.2.sites";

    size_t num_sequences    = 12000;
    size_t seq_length       = 100;
    size_t segment_length   = 1000;
    std::string synthetic_motif = "TATATA";
    double synthetic_injection_rate = 0.45;
    unsigned int seed = 42;

    std::vector<size_t> k_values = {6, 8, 9, 10, 12};
    double significant_z     = 6.0;
    double sub_significant_z = 3.5;
    size_t top_seeds_per_k   = 25;
    size_t max_seed_trials_per_round = 80;

    int max_motifs          = 5;
    int max_em_iterations   = 20;
    int em_patience         = 2;
    double em_min_improvement = 1e-3;
    size_t min_sites_for_em = 8;
    size_t max_sites_for_mstep = 60000;
    double pseudocount = 0.2;
    double min_enrichment_to_report = 1.20;
    double max_mask_fraction = 0.85;
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

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

void print_separator(char ch = '=', int width = 90) {
    std::cout << std::string(width, ch) << '\n';
}

EMConfig make_em_cfg(const PipelineConfig& cfg) {
    EMConfig em;
    em.max_iterations  = cfg.max_em_iterations;
    em.patience        = cfg.em_patience;
    em.min_improvement = cfg.em_min_improvement;
    em.min_sites       = cfg.min_sites_for_em;
    em.max_mstep_sites = cfg.max_sites_for_mstep;
    em.pseudocount     = cfg.pseudocount;
    return em;
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-seed refinement
// ─────────────────────────────────────────────────────────────────────────────

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
    const auto initial_sites = collect_hd_sites(seed.kmer, masked_primary);
    if (initial_sites.size() < cfg.min_sites_for_em) {
        return std::nullopt;
    }

    PWM init_pwm = build_pwm_from_sites(initial_sites, cfg.pseudocount);
    if (init_pwm.width == 0) {
        return std::nullopt;
    }

    const EMConfig em_cfg = make_em_cfg(cfg);
    auto refined = run_em_refinement(
        init_pwm, primary, control,
        primary_mask, control_mask,
        background, em_cfg, seed.kmer
    );

    if (!refined.has_value() ||
        refined->enrichment < cfg.min_enrichment_to_report) {
        return std::nullopt;
    }

    MotifRoundResult result;
    result.seed          = seed.kmer;
    result.consensus     = get_consensus(refined->pwm);
    result.pwm           = refined->pwm;
    result.width         = refined->pwm.width;
    result.z_score       = seed.z_score;
    result.threshold     = refined->learned_threshold;
    result.enrichment    = refined->enrichment;
    result.primary_hits  = refined->primary_hits;
    result.control_hits  = refined->control_hits;
    result.em_iterations = refined->iterations;
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Dataset construction
// ─────────────────────────────────────────────────────────────────────────────

std::pair<SequenceList, SequenceList> build_dataset(const PipelineConfig& cfg) {
    if (cfg.use_real_data) {
        return generate_real_test_data(
            cfg.baseline_file, cfg.sites_file,
            cfg.num_sequences, cfg.segment_length, cfg.seed
        );
    }

    DataGenerator gen(cfg.seed);
    return gen.generate_test_data(
        cfg.num_sequences, cfg.seq_length,
        cfg.synthetic_motif, cfg.synthetic_injection_rate
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Argument parsing
// ─────────────────────────────────────────────────────────────────────────────

PipelineConfig parse_args(int argc, char* argv[]) {
    PipelineConfig cfg;

    auto parse_size_t = [](const std::string& s) {
        size_t pos = 0;
        const size_t v = std::stoull(s, &pos, 10);
        if (pos != s.size()) throw std::invalid_argument("Invalid integer");
        return v;
    };
    auto parse_int = [](const std::string& s) {
        size_t pos = 0;
        const int v = std::stoi(s, &pos, 10);
        if (pos != s.size()) throw std::invalid_argument("Invalid integer");
        return v;
    };
    auto parse_double = [](const std::string& s) {
        size_t pos = 0;
        const double v = std::stod(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("Invalid float");
        return v;
    };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument(std::string("Missing value for ") + name);
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

    if (cfg.synthetic_injection_rate < 0.0) cfg.synthetic_injection_rate = 0.0;
    if (cfg.synthetic_injection_rate > 1.0) cfg.synthetic_injection_rate = 1.0;
    cfg.max_motifs = std::max(1, cfg.max_motifs);
    return cfg;
}

} // namespace

// =============================================================================
// main
// =============================================================================

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

    // Load dataset
    SequenceList primary, control;
    try {
        std::tie(primary, control) = build_dataset(cfg);
    } catch (const std::exception& e) {
        std::cerr << "Dataset error: " << e.what() << '\n';
        return 1;
    }

    std::cout << "Loaded dataset: primary=" << primary.size()
              << ", control=" << control.size() << "\n";
    if (primary.empty() || control.empty()) {
        std::cerr << "Empty dataset; aborting.\n";
        return 1;
    }

    // Train Markov background
    MarkovOrder3 background;
    background.train(control);

    // Initialise masks
    MaskMatrix primary_mask = init_mask(primary);
    MaskMatrix control_mask = init_mask(control);

    // Multi-round discovery with erasing
    std::vector<MotifRoundResult> discovered;
    discovered.reserve(static_cast<size_t>(cfg.max_motifs));

    for (int round = 0; round < cfg.max_motifs; ++round) {
        print_separator('-');
        std::cout << "Round " << (round + 1) << '\n';
        print_separator('-');

        const SequenceList masked_primary = apply_mask(primary, primary_mask);
        const SequenceList masked_control = apply_mask(control, control_mask);

        auto seed_pool = build_seed_pool(
            masked_primary, masked_control,
            cfg.k_values, cfg.sub_significant_z, cfg.top_seeds_per_k
        );
        if (seed_pool.empty()) {
            std::cout << "No candidate seed passed z-thresholds.\n";
            break;
        }

        std::optional<MotifRoundResult> best;
        size_t tried = 0;
        for (const auto& seed : seed_pool) {
            if (tried++ >= cfg.max_seed_trials_per_round) break;

            auto refined = refine_seed_with_em(
                seed, primary, control,
                primary_mask, control_mask,
                background, cfg
            );
            if (!refined.has_value()) continue;

            if (!best.has_value() || refined->enrichment > best->enrichment) {
                best = refined;
            }
        }

        if (!best.has_value()) {
            std::cout << "No seed converged to a valid motif in this round.\n";
            break;
        }

        // Erase discovered motif
        const int erased_primary = erase_by_pwm(
            primary, primary_mask, best->pwm, background, best->threshold
        );
        const int erased_control = erase_by_pwm(
            control, control_mask, best->pwm, background, best->threshold
        );

        std::cout << "Seed                : " << best->seed
                  << " (z=" << std::fixed << std::setprecision(3) << best->z_score << ")\n";
        std::cout << "Consensus           : " << best->consensus
                  << " (W=" << best->width << ")\n";
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
        std::cout << "Primary masked frac : " << std::fixed << std::setprecision(2)
                  << (frac * 100.0) << "%\n";
        if (frac >= cfg.max_mask_fraction) {
            std::cout << "Mask fraction limit reached; stopping.\n";
            break;
        }
    }

    // Summary
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
