#pragma once

#include <cstddef>
#include <cstdint>

#include "helix/physics.hpp"

namespace helix::engine::cuda {

bool is_available();

void score_pairs_encoded(
    const std::uint8_t* guides,
    std::size_t n_guides,
    const std::uint8_t* windows,
    std::size_t n_windows,
    std::size_t seq_len,
    const PhysicsParams& params,
    const float* pam_penalties,
    float* out_scores);

}  // namespace helix::engine::cuda
