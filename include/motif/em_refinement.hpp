/**
 * @file em_refinement.hpp
 * @brief EM-style PWM refinement with Markov background LLR scoring.
 *
 * Provides the iterative E-step / threshold / M-step loop, mask
 * utilities, and a high-level `run_em_refinement` entry point.
 */

#ifndef MOTIF_EM_REFINEMENT_HPP
#define MOTIF_EM_REFINEMENT_HPP

#include "motif/markov.hpp"
#include "motif/types.hpp"

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace motif {

// ─────────────────────────────────────────────────────────────────────────────
// Data structures
// ─────────────────────────────────────────────────────────────────────────────

/** @brief Best-match result for a single sequence. */
struct MatchResult {
    int position = -1;
    double score = -std::numeric_limits<double>::infinity();
};

/** @brief Outcome of threshold optimisation (best enrichment point). */
struct ThresholdResult {
    bool valid = false;
    double threshold = 0.0;
    double enrichment = 0.0;
    int primary_hits = 0;
    int control_hits = 0;
};

/** @brief A fully refined motif model. */
struct RefinedModel {
    PWM pwm;
    double learned_threshold = 0.0;
    double enrichment = -std::numeric_limits<double>::infinity();
    int iterations = 0;
    int primary_hits = 0;
    int control_hits = 0;
    std::string source_seed;
};

/** @brief Configuration knobs for the EM loop. */
struct EMConfig {
    int max_iterations = 20;
    int patience = 2;
    double min_improvement = 1e-3;
    size_t min_sites = 8;
    size_t max_mstep_sites = 60000;
    double pseudocount = 0.2;
};

// ─────────────────────────────────────────────────────────────────────────────
// Mask types and utilities
// ─────────────────────────────────────────────────────────────────────────────

using MaskRow    = std::vector<uint8_t>;
using MaskMatrix = std::vector<MaskRow>;

/** @brief Initialise an all-zero mask matrix matching @p sequences. */
MaskMatrix init_mask(const SequenceList& sequences);

/** @brief Return a copy of @p sequences with masked positions set to 'N'. */
SequenceList apply_mask(const SequenceList& sequences, const MaskMatrix& mask);

/** @brief True if any position in [start, start+width) is masked. */
bool window_is_masked(const MaskRow& mask, size_t start, size_t width);

/** @brief Mark positions [start, start+width) as masked. */
void mark_window(MaskRow& mask, size_t start, size_t width);

/** @brief Fraction of all positions that are masked. */
double masked_fraction(const MaskMatrix& mask);

// ─────────────────────────────────────────────────────────────────────────────
// Scoring functions
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Log-likelihood ratio score for a single window.
 *
 * score = Σ_j log2( PWM[base][j] / bg(seq, start+j) )
 */
double score_window_llr(
    const Sequence& seq,
    size_t start,
    const PWM& pwm,
    const MarkovOrder3& bg
);

/** @brief Best LLR match in a sequence (no mask). */
MatchResult find_best_match_llr(
    const Sequence& seq,
    const PWM& pwm,
    const MarkovOrder3& bg
);

/** @brief Best LLR match in a sequence, skipping masked windows. */
MatchResult find_best_match_llr(
    const Sequence& seq,
    const MaskRow& mask,
    const PWM& pwm,
    const MarkovOrder3& bg
);

/** @brief Scan all sequences (no mask). */
std::vector<MatchResult> scan_best_matches_llr(
    const SequenceList& sequences,
    const PWM& pwm,
    const MarkovOrder3& bg
);

/** @brief Scan all sequences, respecting masks. */
std::vector<MatchResult> scan_best_matches_llr(
    const SequenceList& sequences,
    const MaskMatrix& mask,
    const PWM& pwm,
    const MarkovOrder3& bg
);

// ─────────────────────────────────────────────────────────────────────────────
// Threshold & site collection
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Choose the threshold that maximises enrichment(primary / control).
 */
ThresholdResult choose_threshold(
    const std::vector<MatchResult>& primary_matches,
    const std::vector<MatchResult>& control_matches
);

/**
 * @brief Collect aligned sites from sequences whose match score ≥ threshold.
 */
std::vector<std::string> collect_mstep_sites(
    const SequenceList& sequences,
    const std::vector<MatchResult>& matches,
    double threshold,
    size_t width,
    size_t max_sites
);

// ─────────────────────────────────────────────────────────────────────────────
// EM refinement
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Run the full EM refinement loop starting from @p initial_pwm.
 *
 * The loop iterates: scan → choose_threshold → collect → rebuild PWM,
 * with patience-based early stopping.
 *
 * @param initial_pwm   Starting PWM (e.g. from HD-site collection).
 * @param primary        Primary sequence set.
 * @param control        Control sequence set.
 * @param bg             Markov background model.
 * @param cfg            EM configuration knobs.
 * @param source_seed    Label carried in the returned RefinedModel.
 * @return Best model found, or nullopt if refinement failed.
 */
std::optional<RefinedModel> run_em_refinement(
    const PWM& initial_pwm,
    const SequenceList& primary,
    const SequenceList& control,
    const MarkovOrder3& bg,
    const EMConfig& cfg,
    const std::string& source_seed = ""
);

/**
 * @brief Overload that accepts mask matrices (used by the erasing loop).
 */
std::optional<RefinedModel> run_em_refinement(
    const PWM& initial_pwm,
    const SequenceList& primary,
    const SequenceList& control,
    const MaskMatrix& primary_mask,
    const MaskMatrix& control_mask,
    const MarkovOrder3& bg,
    const EMConfig& cfg,
    const std::string& source_seed = ""
);

// ─────────────────────────────────────────────────────────────────────────────
// Erasing helper
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Mask all windows scoring ≥ threshold against a PWM.
 * @return Number of newly masked windows.
 */
int erase_by_pwm(
    const SequenceList& sequences,
    MaskMatrix& mask,
    const PWM& pwm,
    const MarkovOrder3& bg,
    double threshold
);

} // namespace motif

#endif // MOTIF_EM_REFINEMENT_HPP
