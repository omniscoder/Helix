#pragma once

#include <cstddef>
#include <cstdint>

namespace helix::engine {

struct PhysicsParams {
    const float* weights;
    std::size_t weights_len;
    float bulge_penalty;
    float gc_low;
    float gc_high;
    float pam_penalty;
    float min_score;
};

// Score a single (guide, window) pair. Arrays use the same uint8 encoding as
// helix.crispr.physics._encode_sequence_to_uint8.
float compute_on_target_score_encoded(
    const std::uint8_t* guide,
    std::size_t guide_len,
    const std::uint8_t* window,
    std::size_t window_len,
    const PhysicsParams& params);

// Batched kernel contract: guides is (G, L) row-major, windows is (N, L)
// row-major, and out_scores is (G, N) row-major. Shapes must match exactly.
void score_pairs_encoded(
    const std::uint8_t* guides,
    std::size_t n_guides,
    const std::uint8_t* windows,
    std::size_t n_windows,
    std::size_t seq_len,
    const PhysicsParams& params,
    const float* pam_penalties,
    float* out_scores);

}  // namespace helix::engine
