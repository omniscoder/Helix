#include "helix/physics.hpp"

#include <algorithm>
#include <cmath>

namespace helix::engine {

namespace {

inline bool is_gc(std::uint8_t base) {
    return base == 1 || base == 2;  // C or G
}

float gc_penalty(const std::uint8_t* window, std::size_t len, float gc_low, float gc_high) {
    if (len == 0) {
        return 0.0f;
    }
    std::size_t gc_count = 0;
    for (std::size_t idx = 0; idx < len; ++idx) {
        if (is_gc(window[idx])) {
            ++gc_count;
        }
    }
    const float gc_fraction = static_cast<float>(gc_count) / static_cast<float>(len);
    if (gc_fraction >= gc_low && gc_fraction <= gc_high) {
        return 0.0f;
    }
    const float delta = std::min(std::abs(gc_fraction - gc_low), std::abs(gc_fraction - gc_high));
    return delta * 2.0f;
}

float mismatch_cost(
    const std::uint8_t* guide,
    std::size_t guide_len,
    const std::uint8_t* window,
    std::size_t window_len,
    const PhysicsParams& params) {
    const std::size_t overlap = std::min(guide_len, window_len);
    float cost = 0.0f;
    for (std::size_t idx = 0; idx < overlap; ++idx) {
        if (guide[idx] != window[idx]) {
            const float weight = idx < params.weights_len ? params.weights[idx] : 1.0f;
            cost += weight;
        }
    }
    if (window_len > guide_len) {
        const float extra = static_cast<float>(window_len - guide_len);
        cost += 2.0f * extra * params.bulge_penalty;
    } else if (guide_len > window_len) {
        const float extra = static_cast<float>(guide_len - window_len);
        cost += extra * params.bulge_penalty;
    }
    return cost;
}

float compute_score(
    float mismatch,
    float gc_pen,
    float pam_pen,
    float guide_len,
    float min_score) {
    const float span = std::max(1.0f, guide_len);
    const float total = mismatch + gc_pen + pam_pen;
    const float raw = std::max(min_score, 1.0f - (total / span));
    return std::min(1.0f, raw);
}

}  // namespace

float compute_on_target_score_encoded(
    const std::uint8_t* guide,
    std::size_t guide_len,
    const std::uint8_t* window,
    std::size_t window_len,
    const PhysicsParams& params) {
    const float mismatch = mismatch_cost(guide, guide_len, window, window_len, params);
    const float gc = gc_penalty(window, window_len, params.gc_low, params.gc_high);
    return compute_score(mismatch, gc, params.pam_penalty, static_cast<float>(guide_len), params.min_score);
}

void score_pairs_encoded(
    const std::uint8_t* guides,
    std::size_t n_guides,
    const std::uint8_t* windows,
    std::size_t n_windows,
    std::size_t seq_len,
    const PhysicsParams& params,
    const float* pam_penalties,
    float* out_scores) {
    for (std::size_t g = 0; g < n_guides; ++g) {
        const std::uint8_t* guide = guides + (g * seq_len);
        for (std::size_t w = 0; w < n_windows; ++w) {
            const std::uint8_t* window = windows + (w * seq_len);
            PhysicsParams local = params;
            if (pam_penalties != nullptr) {
                local.pam_penalty = pam_penalties[g * n_windows + w];
            }
            out_scores[g * n_windows + w] = compute_on_target_score_encoded(
                guide,
                seq_len,
                window,
                seq_len,
                local);
        }
    }
}

}  // namespace helix::engine
