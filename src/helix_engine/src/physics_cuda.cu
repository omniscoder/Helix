#include "helix/physics_cuda.hpp"

#ifdef HELIX_ENGINE_HAS_CUDA

#include <cuda_runtime.h>

#include <cmath>
#include <stdexcept>

namespace helix::engine::cuda {

namespace {

inline void check_cuda(cudaError_t err, const char* context) {
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA error during ") + context + ": " + cudaGetErrorString(err));
    }
}

__device__ float mismatch_cost(
    const std::uint8_t* guide,
    const std::uint8_t* window,
    std::size_t seq_len,
    const float* weights,
    std::size_t weights_len,
    float bulge_penalty) {
    const std::size_t overlap = seq_len;
    float cost = 0.0f;
    for (std::size_t idx = 0; idx < overlap; ++idx) {
        if (guide[idx] != window[idx]) {
            float weight = 1.0f;
            if (idx < weights_len) {
                weight = weights[idx];
            }
            cost += weight;
        }
    }
    // sequences are the same length (seq_len), so no bulge terms beyond mismatch weights
    return cost;
}

__device__ float gc_penalty(const std::uint8_t* window, std::size_t len, float gc_low, float gc_high) {
    if (len == 0) {
        return 0.0f;
    }
    std::size_t gc_count = 0;
    for (std::size_t idx = 0; idx < len; ++idx) {
        if (window[idx] == 1 || window[idx] == 2) {  // C or G
            ++gc_count;
        }
    }
    const float gc_frac = static_cast<float>(gc_count) / static_cast<float>(len);
    if (gc_frac >= gc_low && gc_frac <= gc_high) {
        return 0.0f;
    }
    const float delta = fminf(fabsf(gc_frac - gc_low), fabsf(gc_frac - gc_high));
    return delta * 2.0f;
}

__device__ float compute_score(float mismatch, float gc_pen, float pam_pen, float guide_len, float min_score) {
    const float span = fmaxf(1.0f, guide_len);
    const float total = mismatch + gc_pen + pam_pen;
    const float raw = fmaxf(min_score, 1.0f - (total / span));
    return fminf(1.0f, raw);
}

__global__ void score_pairs_kernel(
    const std::uint8_t* guides,
    std::size_t n_guides,
    const std::uint8_t* windows,
    std::size_t n_windows,
    std::size_t seq_len,
    PhysicsParams params,
    const float* pam_penalties,
    float* out_scores) {
    const std::size_t w = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t g = blockIdx.y * blockDim.y + threadIdx.y;
    if (g >= n_guides || w >= n_windows) {
        return;
    }
    const std::uint8_t* guide = guides + g * seq_len;
    const std::uint8_t* window = windows + w * seq_len;
    const float pam_pen = pam_penalties ? pam_penalties[g * n_windows + w] : 0.0f;
    const float mismatch = mismatch_cost(guide, window, seq_len, params.weights, params.weights_len, params.bulge_penalty);
    const float gc = gc_penalty(window, seq_len, params.gc_low, params.gc_high);
    out_scores[g * n_windows + w] = compute_score(mismatch, gc, pam_pen, static_cast<float>(seq_len), params.min_score);
}

}  // namespace

bool is_available() {
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess) {
        return false;
    }
    return count > 0;
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
    if (!is_available()) {
        throw std::runtime_error("No CUDA devices available for Helix engine");
    }
    const std::size_t guides_bytes = n_guides * seq_len * sizeof(std::uint8_t);
    const std::size_t windows_bytes = n_windows * seq_len * sizeof(std::uint8_t);
    const std::size_t scores_bytes = n_guides * n_windows * sizeof(float);
    const std::size_t weights_bytes = params.weights_len * sizeof(float);
    const std::size_t pam_bytes = n_guides * n_windows * sizeof(float);

    std::uint8_t* d_guides = nullptr;
    std::uint8_t* d_windows = nullptr;
    float* d_scores = nullptr;
    float* d_weights = nullptr;
    float* d_pams = nullptr;

    check_cuda(cudaMalloc(&d_guides, guides_bytes), "cudaMalloc guides");
    check_cuda(cudaMalloc(&d_windows, windows_bytes), "cudaMalloc windows");
    check_cuda(cudaMalloc(&d_scores, scores_bytes), "cudaMalloc scores");
    check_cuda(cudaMalloc(&d_weights, weights_bytes), "cudaMalloc weights");
    if (pam_penalties) {
        check_cuda(cudaMalloc(&d_pams, pam_bytes), "cudaMalloc pam");
    }

    check_cuda(cudaMemcpy(d_guides, guides, guides_bytes, cudaMemcpyHostToDevice), "memcpy guides");
    check_cuda(cudaMemcpy(d_windows, windows, windows_bytes, cudaMemcpyHostToDevice), "memcpy windows");
    check_cuda(cudaMemcpy(d_weights, params.weights, weights_bytes, cudaMemcpyHostToDevice), "memcpy weights");
    if (pam_penalties && d_pams) {
        check_cuda(cudaMemcpy(d_pams, pam_penalties, pam_bytes, cudaMemcpyHostToDevice), "memcpy pam");
    }

    PhysicsParams device_params = params;
    device_params.weights = d_weights;

    dim3 block(16, 16);
    dim3 grid((n_windows + block.x - 1) / block.x, (n_guides + block.y - 1) / block.y);
    score_pairs_kernel<<<grid, block>>>(
        d_guides,
        n_guides,
        d_windows,
        n_windows,
        seq_len,
        device_params,
        d_pams,
        d_scores);
    check_cuda(cudaGetLastError(), "score_pairs_kernel launch");
    check_cuda(cudaDeviceSynchronize(), "score_pairs_kernel sync");

    check_cuda(cudaMemcpy(out_scores, d_scores, scores_bytes, cudaMemcpyDeviceToHost), "memcpy scores");

    cudaFree(d_guides);
    cudaFree(d_windows);
    cudaFree(d_scores);
    cudaFree(d_weights);
    if (d_pams) {
        cudaFree(d_pams);
    }
}

}  // namespace helix::engine::cuda

#else

#include <stdexcept>

namespace helix::engine::cuda {

bool is_available() {
    return false;
}

void score_pairs_encoded(
    const std::uint8_t*,
    std::size_t,
    const std::uint8_t*,
    std::size_t,
    std::size_t,
    const PhysicsParams&,
    const float*,
    float*) {
    throw std::runtime_error("Helix engine built without CUDA support");
}

}  // namespace helix::engine::cuda

#endif
