/**
 * @file em_refinement.cpp
 * @brief Implementation of EM-style PWM refinement with Markov LLR scoring.
 */

#include "motif/em_refinement.hpp"
#include "motif/pwm_construction.hpp"

#include <algorithm>
#include <cmath>

namespace motif {

// ─────────────────────────────────────────────────────────────────────────────
// Mask utilities
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// Scoring
// ─────────────────────────────────────────────────────────────────────────────

double score_window_llr(
    const Sequence& seq,
    size_t start,
    const PWM& pwm,
    const MarkovOrder3& bg
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
        const double bg_prob  = std::max(bg.conditional_prob(seq, pos), kMinProb);
        score += std::log2(pwm_prob / bg_prob);
    }

    return score;
}

MatchResult find_best_match_llr(
    const Sequence& seq,
    const PWM& pwm,
    const MarkovOrder3& bg
) {
    MatchResult best;
    if (pwm.width == 0 || seq.size() < pwm.width) {
        return best;
    }

    for (size_t start = 0; start + pwm.width <= seq.size(); ++start) {
        const double s = score_window_llr(seq, start, pwm, bg);
        if (s > best.score) {
            best.score = s;
            best.position = static_cast<int>(start);
        }
    }

    return best;
}

MatchResult find_best_match_llr(
    const Sequence& seq,
    const MaskRow& mask,
    const PWM& pwm,
    const MarkovOrder3& bg
) {
    MatchResult best;
    if (pwm.width == 0 || seq.size() < pwm.width) {
        return best;
    }

    for (size_t start = 0; start + pwm.width <= seq.size(); ++start) {
        if (window_is_masked(mask, start, pwm.width)) {
            continue;
        }
        const double s = score_window_llr(seq, start, pwm, bg);
        if (s > best.score) {
            best.score = s;
            best.position = static_cast<int>(start);
        }
    }

    return best;
}

std::vector<MatchResult> scan_best_matches_llr(
    const SequenceList& sequences,
    const PWM& pwm,
    const MarkovOrder3& bg
) {
    std::vector<MatchResult> results;
    results.reserve(sequences.size());

    for (const auto& seq : sequences) {
        results.push_back(find_best_match_llr(seq, pwm, bg));
    }

    return results;
}

std::vector<MatchResult> scan_best_matches_llr(
    const SequenceList& sequences,
    const MaskMatrix& mask,
    const PWM& pwm,
    const MarkovOrder3& bg
) {
    std::vector<MatchResult> results;
    results.reserve(sequences.size());

    const size_t n = std::min(sequences.size(), mask.size());
    for (size_t i = 0; i < n; ++i) {
        results.push_back(find_best_match_llr(sequences[i], mask[i], pwm, bg));
    }

    return results;
}

// ─────────────────────────────────────────────────────────────────────────────
// Threshold & site collection
// ─────────────────────────────────────────────────────────────────────────────

ThresholdResult choose_threshold(
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
    candidates.erase(
        std::unique(candidates.begin(), candidates.end()),
        candidates.end()
    );

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

        const double p_rate =
            (static_cast<double>(hp) + 1.0) /
            (static_cast<double>(primary_matches.size()) + 2.0);
        const double c_rate =
            (static_cast<double>(hc) + 1.0) /
            (static_cast<double>(control_matches.size()) + 2.0);
        const double enrichment = p_rate / std::max(c_rate, kMinProb);

        if (!best.valid || enrichment > best.enrichment ||
            (std::abs(enrichment - best.enrichment) < 1e-9 &&
             hp > best.primary_hits)) {
            best.valid = true;
            best.threshold = threshold;
            best.enrichment = enrichment;
            best.primary_hits = hp;
            best.control_hits = hc;
        }
    }

    return best;
}

std::vector<std::string> collect_mstep_sites(
    const SequenceList& sequences,
    const std::vector<MatchResult>& matches,
    double threshold,
    size_t width,
    size_t max_sites
) {
    std::vector<std::string> sites;
    sites.reserve(std::min(sequences.size(), max_sites));

    const size_t n = std::min(sequences.size(), matches.size());
    for (size_t i = 0; i < n && sites.size() < max_sites; ++i) {
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

// ─────────────────────────────────────────────────────────────────────────────
// EM refinement – no mask
// ─────────────────────────────────────────────────────────────────────────────

std::optional<RefinedModel> run_em_refinement(
    const PWM& initial_pwm,
    const SequenceList& primary,
    const SequenceList& control,
    const MarkovOrder3& bg,
    const EMConfig& cfg,
    const std::string& source_seed
) {
    if (initial_pwm.width == 0) {
        return std::nullopt;
    }

    PWM pwm = initial_pwm;
    RefinedModel best;
    best.source_seed = source_seed;

    double prev_enrichment = -std::numeric_limits<double>::infinity();
    int stale = 0;

    for (int iter = 0; iter < cfg.max_iterations; ++iter) {
        const auto pm = scan_best_matches_llr(primary, pwm, bg);
        const auto cm = scan_best_matches_llr(control, pwm, bg);

        const ThresholdResult thr = choose_threshold(pm, cm);
        if (!thr.valid ||
            thr.primary_hits < static_cast<int>(cfg.min_sites)) {
            break;
        }

        const auto sites = collect_mstep_sites(
            primary, pm, thr.threshold, pwm.width, cfg.max_mstep_sites
        );
        if (sites.size() < cfg.min_sites) {
            break;
        }

        const PWM next_pwm = build_pwm_from_sites(sites, cfg.pseudocount);
        if (next_pwm.width == 0) {
            break;
        }

        if (thr.enrichment > best.enrichment) {
            best.pwm = next_pwm;
            best.learned_threshold = thr.threshold;
            best.enrichment = thr.enrichment;
            best.iterations = iter + 1;
            best.primary_hits = thr.primary_hits;
            best.control_hits = thr.control_hits;
        }

        if (prev_enrichment > -std::numeric_limits<double>::infinity() &&
            thr.enrichment - prev_enrichment < cfg.min_improvement) {
            ++stale;
        } else {
            stale = 0;
        }

        prev_enrichment = thr.enrichment;
        pwm = next_pwm;

        if (stale >= cfg.patience) {
            break;
        }
    }

    if (best.pwm.width == 0) {
        return std::nullopt;
    }
    return best;
}

// ─────────────────────────────────────────────────────────────────────────────
// EM refinement – with mask (for erasing loop)
// ─────────────────────────────────────────────────────────────────────────────

std::optional<RefinedModel> run_em_refinement(
    const PWM& initial_pwm,
    const SequenceList& primary,
    const SequenceList& control,
    const MaskMatrix& primary_mask,
    const MaskMatrix& control_mask,
    const MarkovOrder3& bg,
    const EMConfig& cfg,
    const std::string& source_seed
) {
    if (initial_pwm.width == 0) {
        return std::nullopt;
    }

    PWM pwm = initial_pwm;
    RefinedModel best;
    best.source_seed = source_seed;

    double prev_enrichment = -std::numeric_limits<double>::infinity();
    int stale = 0;

    for (int iter = 0; iter < cfg.max_iterations; ++iter) {
        const auto pm = scan_best_matches_llr(primary, primary_mask, pwm, bg);
        const auto cm = scan_best_matches_llr(control, control_mask, pwm, bg);

        const ThresholdResult thr = choose_threshold(pm, cm);
        if (!thr.valid ||
            thr.primary_hits < static_cast<int>(cfg.min_sites)) {
            break;
        }

        const auto sites = collect_mstep_sites(
            primary, pm, thr.threshold, pwm.width, cfg.max_mstep_sites
        );
        if (sites.size() < cfg.min_sites) {
            break;
        }

        const PWM next_pwm = build_pwm_from_sites(sites, cfg.pseudocount);
        if (next_pwm.width == 0) {
            break;
        }

        if (thr.enrichment > best.enrichment) {
            best.pwm = next_pwm;
            best.learned_threshold = thr.threshold;
            best.enrichment = thr.enrichment;
            best.iterations = iter + 1;
            best.primary_hits = thr.primary_hits;
            best.control_hits = thr.control_hits;
        }

        if (prev_enrichment > -std::numeric_limits<double>::infinity() &&
            thr.enrichment - prev_enrichment < cfg.min_improvement) {
            ++stale;
        } else {
            stale = 0;
        }

        prev_enrichment = thr.enrichment;
        pwm = next_pwm;

        if (stale >= cfg.patience) {
            break;
        }
    }

    if (best.pwm.width == 0) {
        return std::nullopt;
    }
    return best;
}

// ─────────────────────────────────────────────────────────────────────────────
// Erasing helper
// ─────────────────────────────────────────────────────────────────────────────

int erase_by_pwm(
    const SequenceList& sequences,
    MaskMatrix& mask,
    const PWM& pwm,
    const MarkovOrder3& bg,
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
            const double score = score_window_llr(seq, start, pwm, bg);
            if (score >= threshold) {
                mark_window(mask[i], start, pwm.width);
                ++hit_count;
            }
        }
    }

    return hit_count;
}

} // namespace motif
