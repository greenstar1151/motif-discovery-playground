/**
 * @file morbius_eval.cpp
 * @brief Morbius benchmark evaluator — thin CLI wrapper around motif_core.
 *
 * Pipeline:
 *   1) Load input FASTA + build/load control set
 *   2) Multi-k seed discovery with two-proportion z-test
 *   3) HD≤1 site collection for initialisation
 *   4) EM-style PWM refinement with Markov order-3 background LLR scoring
 *   5) Best-site prediction per sequence → Morbius TSV output
 */

#include "motif.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace motif;
using namespace std::chrono;

namespace {

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

struct Config {
    std::string input_file;
    std::string output_file;
    std::string control_file;
    std::string seed_motif;
    std::string jaspar_file;

    size_t k = 16;
    size_t motif_width = 0;
    size_t top_seeds = 8;

    double threshold = 0.0;
    unsigned int rng_seed = 42;

    double z_sig = 6.0;
    double z_sub = 3.5;
    int em_iters = 20;
    int em_patience = 2;
    double em_min_improve = 1e-3;
    size_t min_sites = 8;
    size_t max_mstep_sites = 60000;
    double pseudocount = 0.2;

    std::string multi_k_csv;
};

struct MotifHit {
    size_t seq_id = 0;
    int position = -1;
    std::string motif;
    double score = -std::numeric_limits<double>::infinity();
};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

void print_usage(const char* prog) {
    std::cerr
        << "\nUsage: " << prog << " <input.fasta> <output.tsv> [options]\n\n"
        << "Options:\n"
        << "  -k <int>              Base seed length (default: 16)\n"
        << "  -m <int>              Motif width (default: same as k)\n"
        << "  --control <file>      Control FASTA file (default: shuffled input)\n"
        << "  --seed <string>       Use a fixed seed motif\n"
        << "  --pwm <file>          Load PWM from JASPAR file\n"
        << "  --top <int>           Number of seed candidates to refine (default: 8)\n"
        << "  --threshold <float>   Output threshold (default: 0.0 => always best site)\n"
        << "  --multi-k <list>      Seed widths (e.g. 10,12,14,16)\n"
        << "  --zsig <float>        Significant z cutoff (default: 6.0)\n"
        << "  --zsub <float>        Sub-significant z cutoff (default: 3.5)\n"
        << "  --em-iters <int>      Max EM refinement iterations (default: 20)\n"
        << "  --rng-seed <int>      RNG seed (default: 42)\n"
        << "  --help                Show this help\n\n"
        << "Output format (TSV):\n"
        << "  seq_id    position    motif\n\n"
        << "Example:\n"
        << "  " << prog
        << " DATASET_DNA_3.fasta DATASET_DNA_3_result.tsv -k 16 --top 12\n";
}

EMConfig make_em_cfg(const Config& cfg) {
    EMConfig em;
    em.max_iterations  = cfg.em_iters;
    em.patience        = cfg.em_patience;
    em.min_improvement = cfg.em_min_improve;
    em.min_sites       = cfg.min_sites;
    em.max_mstep_sites = cfg.max_mstep_sites;
    em.pseudocount     = cfg.pseudocount;
    return em;
}

std::vector<size_t> derive_k_values(const Config& cfg) {
    if (!cfg.multi_k_csv.empty()) {
        std::set<size_t> custom;
        std::stringstream ss(cfg.multi_k_csv);
        std::string token;

        while (std::getline(ss, token, ',')) {
            if (token.empty()) continue;
            const size_t v = static_cast<size_t>(std::stoull(token));
            if (v >= 4 && v <= MAX_MOTIF_WIDTH) {
                custom.insert(v);
            }
        }
        if (!custom.empty()) {
            return std::vector<size_t>(custom.begin(), custom.end());
        }
    }

    const size_t base = std::max<size_t>(4, std::min(cfg.k, MAX_MOTIF_WIDTH));
    std::set<size_t> values;
    values.insert(base);
    if (base > 6)  values.insert(base - 2);
    if (base > 8)  values.insert(base - 4);
    if (base + 2 <= MAX_MOTIF_WIDTH) values.insert(base + 2);
    if (cfg.motif_width > 0 && cfg.motif_width <= MAX_MOTIF_WIDTH) {
        values.insert(cfg.motif_width);
    }

    return std::vector<size_t>(values.begin(), values.end());
}

// ─────────────────────────────────────────────────────────────────────────────
// Single-seed refinement
// ─────────────────────────────────────────────────────────────────────────────

std::optional<RefinedModel> refine_from_seed(
    const std::string& seed,
    const SequenceList& primary,
    const SequenceList& control,
    const MarkovOrder3& bg,
    const Config& cfg
) {
    const size_t target_width = (cfg.motif_width > 0) ? cfg.motif_width : seed.size();
    const auto init_sites = collect_hd_sites(seed, primary, target_width, 1);
    if (init_sites.size() < cfg.min_sites) {
        return std::nullopt;
    }

    const PWM init_pwm = build_pwm_from_sites(init_sites, cfg.pseudocount);
    if (init_pwm.width == 0) {
        return std::nullopt;
    }

    return run_em_refinement(init_pwm, primary, control, bg, make_em_cfg(cfg), seed);
}

// ─────────────────────────────────────────────────────────────────────────────
// Auto-discovery: build seed pool → refine top seeds → pick best
// ─────────────────────────────────────────────────────────────────────────────

std::optional<RefinedModel> discover_and_refine(
    const SequenceList& primary,
    const SequenceList& control,
    const MarkovOrder3& bg,
    const Config& cfg
) {
    const auto k_values = derive_k_values(cfg);
    const auto pool = build_seed_pool(primary, control, k_values, cfg.z_sub);
    if (pool.empty()) {
        return std::nullopt;
    }

    std::optional<RefinedModel> best;
    const size_t trials = std::min(pool.size(), std::max<size_t>(1, cfg.top_seeds));

    for (size_t i = 0; i < trials; ++i) {
        auto refined = refine_from_seed(pool[i].kmer, primary, control, bg, cfg);
        if (!refined.has_value()) continue;

        if (!best.has_value() || refined->enrichment > best->enrichment) {
            best = refined;
        }
    }

    // Fallback: use top seed with simple PWM if EM failed for all
    if (!best.has_value()) {
        const std::string fallback_seed = pool.front().kmer;
        const PWM fallback_pwm = build_pwm_from_seed(fallback_seed, primary, cfg.pseudocount);
        if (fallback_pwm.width == 0) return std::nullopt;

        RefinedModel fallback;
        fallback.pwm = fallback_pwm;
        fallback.learned_threshold = 0.0;
        fallback.enrichment = 1.0;
        fallback.source_seed = fallback_seed;
        return fallback;
    }

    return best;
}

// ─────────────────────────────────────────────────────────────────────────────
// Discovery orchestrator (3 modes: JASPAR / fixed seed / auto)
// ─────────────────────────────────────────────────────────────────────────────

std::vector<MotifHit> run_discovery(
    const FastaData& input_data,
    const SequenceList& control_seqs,
    const Config& cfg
) {
    std::vector<MotifHit> results;
    results.reserve(input_data.sequences.size());

    MarkovOrder3 bg;
    bg.train(control_seqs);

    std::optional<RefinedModel> model;

    // ── Mode 1: JASPAR PWM ──────────────────────────────────────────────
    if (!cfg.jaspar_file.empty()) {
        std::cerr << "[Info] Loading PWM from JASPAR: " << cfg.jaspar_file << "\n";
        const PWM loaded = load_jaspar_pwm(cfg.jaspar_file);
        model = run_em_refinement(
            loaded, input_data.sequences, control_seqs, bg,
            make_em_cfg(cfg), "JASPAR"
        );
        if (!model.has_value()) {
            RefinedModel fallback;
            fallback.pwm = loaded;
            fallback.source_seed = get_consensus(loaded);
            fallback.enrichment = 1.0;
            model = fallback;
        }
    }
    // ── Mode 2: Fixed seed ──────────────────────────────────────────────
    else if (!cfg.seed_motif.empty()) {
        std::cerr << "[Info] Using fixed seed: " << cfg.seed_motif << "\n";
        model = refine_from_seed(
            cfg.seed_motif, input_data.sequences, control_seqs, bg, cfg
        );
        if (!model.has_value()) {
            const PWM fallback_pwm = build_pwm_from_seed(
                cfg.seed_motif, input_data.sequences, cfg.pseudocount
            );
            if (fallback_pwm.width > 0) {
                RefinedModel fallback;
                fallback.pwm = fallback_pwm;
                fallback.source_seed = cfg.seed_motif;
                fallback.enrichment = 1.0;
                model = fallback;
            }
        }
    }
    // ── Mode 3: Auto-discovery ──────────────────────────────────────────
    else {
        std::cerr << "[Info] Discovering seeds with z-test (base k=" << cfg.k << ")\n";
        model = discover_and_refine(
            input_data.sequences, control_seqs, bg, cfg
        );
    }

    if (!model.has_value() || model->pwm.width == 0) {
        std::cerr << "[Error] Failed to build a valid PWM model\n";
        return results;
    }

    // Report model
    std::cerr << "[Info] Final seed/model: " << model->source_seed << "\n";
    std::cerr << "[Info] Final consensus : " << get_consensus(model->pwm)
              << " (W=" << model->pwm.width << ")\n";
    std::cerr << "[Info] Learned enrich. : " << std::fixed << std::setprecision(3)
              << model->enrichment << "\n";
    std::cerr << "[Info] Learned thres.  : " << std::fixed << std::setprecision(3)
              << model->learned_threshold << "\n";

    // Scan all sequences
    const bool force_threshold = (cfg.threshold != 0.0);
    const double out_threshold = force_threshold ? cfg.threshold : model->learned_threshold;

    std::cerr << "[Info] Scanning " << input_data.sequences.size() << " sequences...\n";
    for (size_t i = 0; i < input_data.sequences.size(); ++i) {
        const auto& seq = input_data.sequences[i];
        const MatchResult best = find_best_match_llr(seq, model->pwm, bg);

        MotifHit hit;
        hit.seq_id = i;
        hit.score  = best.score;

        if (best.position >= 0 && (!force_threshold || best.score >= out_threshold)) {
            hit.position = best.position;
            const size_t pos = static_cast<size_t>(best.position);
            if (pos + model->pwm.width <= seq.size()) {
                hit.motif = seq.substr(pos, model->pwm.width);
            } else if (pos < seq.size()) {
                hit.motif = seq.substr(pos);
            }
        }

        results.push_back(std::move(hit));

        if ((i + 1) % 10000 == 0) {
            std::cerr << "[Progress] " << (i + 1) << " / "
                      << input_data.sequences.size() << "\n";
        }
    }

    return results;
}

// ─────────────────────────────────────────────────────────────────────────────
// TSV writer
// ─────────────────────────────────────────────────────────────────────────────

void write_results_tsv(
    const std::string& filepath,
    const std::vector<MotifHit>& results,
    const std::vector<std::string>& headers
) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + filepath);
    }

    out << "seq_id\tposition\tmotif\n";

    for (const auto& hit : results) {
        std::string seq_id;
        if (hit.seq_id < headers.size()) {
            const std::string& h = headers[hit.seq_id];
            const size_t p = h.find_first_of(" \t");
            seq_id = (p == std::string::npos) ? h : h.substr(0, p);
        } else {
            seq_id = std::to_string(hit.seq_id);
        }

        if (hit.position >= 0) {
            out << seq_id << "\t" << hit.position << "\t" << hit.motif << "\n";
        } else {
            out << seq_id << "\t-1\t\n";
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Argument parsing
// ─────────────────────────────────────────────────────────────────────────────

Config parse_args(int argc, char* argv[]) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        }
    }

    if (argc < 3) {
        print_usage(argv[0]);
        throw std::invalid_argument("Input/output are required");
    }

    Config cfg;
    int i = 1;

    while (i < argc) {
        const std::string arg = argv[i];

        if (arg == "-k" && i + 1 < argc) {
            cfg.k = static_cast<size_t>(std::stoull(argv[++i]));
        } else if (arg == "-m" && i + 1 < argc) {
            cfg.motif_width = static_cast<size_t>(std::stoull(argv[++i]));
        } else if (arg == "--control" && i + 1 < argc) {
            cfg.control_file = argv[++i];
        } else if (arg == "--seed" && i + 1 < argc) {
            cfg.seed_motif = argv[++i];
        } else if (arg == "--pwm" && i + 1 < argc) {
            cfg.jaspar_file = argv[++i];
        } else if (arg == "--top" && i + 1 < argc) {
            cfg.top_seeds = static_cast<size_t>(std::stoull(argv[++i]));
        } else if (arg == "--threshold" && i + 1 < argc) {
            cfg.threshold = std::stod(argv[++i]);
        } else if (arg == "--rng-seed" && i + 1 < argc) {
            cfg.rng_seed = static_cast<unsigned int>(std::stoi(argv[++i]));
        } else if (arg == "--multi-k" && i + 1 < argc) {
            cfg.multi_k_csv = argv[++i];
        } else if (arg == "--zsig" && i + 1 < argc) {
            cfg.z_sig = std::stod(argv[++i]);
        } else if (arg == "--zsub" && i + 1 < argc) {
            cfg.z_sub = std::stod(argv[++i]);
        } else if (arg == "--em-iters" && i + 1 < argc) {
            cfg.em_iters = std::stoi(argv[++i]);
        } else if (cfg.input_file.empty()) {
            cfg.input_file = arg;
        } else if (cfg.output_file.empty()) {
            cfg.output_file = arg;
        } else {
            throw std::invalid_argument("Unknown argument: " + arg);
        }

        ++i;
    }

    if (cfg.input_file.empty() || cfg.output_file.empty()) {
        throw std::invalid_argument("Input and output files are required");
    }

    if (cfg.top_seeds == 0) cfg.top_seeds = 1;
    if (cfg.em_iters <= 0) cfg.em_iters = 1;
    if (cfg.z_sub > cfg.z_sig) cfg.z_sub = cfg.z_sig;
    if (cfg.motif_width > MAX_MOTIF_WIDTH) cfg.motif_width = MAX_MOTIF_WIDTH;
    if (cfg.k > MAX_MOTIF_WIDTH) cfg.k = MAX_MOTIF_WIDTH;

    return cfg;
}

} // namespace

// =============================================================================
// main
// =============================================================================

int main(int argc, char* argv[]) {
    try {
        const Config cfg = parse_args(argc, argv);
        const auto t0 = high_resolution_clock::now();

        std::cerr << "=========================================================\n";
        std::cerr << "  Morbius Benchmark Evaluation (Integrated Pipeline)\n";
        std::cerr << "=========================================================\n\n";

        std::cerr << "[Info] Loading input FASTA: " << cfg.input_file << "\n";
        const FastaData input_data = read_fasta(cfg.input_file);
        if (input_data.sequences.empty()) {
            throw std::runtime_error("No sequences found in input FASTA");
        }

        std::cerr << "[Info] Sequences loaded : " << input_data.sequences.size() << "\n";
        std::cerr << "[Info] Avg length       : " << std::fixed << std::setprecision(1)
                  << input_data.avg_length() << " bp\n";

        SequenceList control;
        if (!cfg.control_file.empty()) {
            std::cerr << "[Info] Loading control FASTA: " << cfg.control_file << "\n";
            FastaData control_data = read_fasta(cfg.control_file);
            control = std::move(control_data.sequences);
        } else {
            std::cerr << "[Info] Generating shuffled control (k-mer preserving, k=2)\n";
            control = generate_control_sequences(input_data.sequences, cfg.rng_seed, 2);
        }

        if (control.empty()) {
            throw std::runtime_error("Control set is empty");
        }
        std::cerr << "[Info] Control sequences: " << control.size() << "\n\n";

        const auto results = run_discovery(input_data, control, cfg);
        if (results.empty()) {
            throw std::runtime_error("No prediction results generated");
        }

        std::cerr << "\n[Info] Writing TSV: " << cfg.output_file << "\n";
        write_results_tsv(cfg.output_file, results, input_data.headers);

        int hits = 0;
        for (const auto& r : results) {
            if (r.position >= 0) ++hits;
        }

        const auto t1 = high_resolution_clock::now();
        const double sec = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

        std::cerr << "\n=========================================================\n";
        std::cerr << "  Summary\n";
        std::cerr << "=========================================================\n";
        std::cerr << "  Total sequences : " << results.size() << "\n";
        std::cerr << "  Motifs found    : " << hits << "\n";
        std::cerr << "  Coverage        : " << std::fixed << std::setprecision(2)
                  << (100.0 * static_cast<double>(hits) /
                      std::max<size_t>(1, results.size()))
                  << "%\n";
        std::cerr << "  Elapsed time    : " << std::fixed << std::setprecision(2)
                  << sec << " sec\n";
        std::cerr << "  Output file     : " << cfg.output_file << "\n\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[Error] " << e.what() << "\n";
        return 1;
    }
}
