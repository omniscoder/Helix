#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "helix/physics.hpp"
#include "helix/physics_cuda.hpp"

namespace py = pybind11;

namespace helix::engine {

namespace {

template <typename T>
py::buffer_info ensure_contiguous(py::array_t<T>& array, const char* name, int expected_dim) {
    if (array.ndim() != expected_dim) {
        throw py::value_error(std::string{name} + " must be " + std::to_string(expected_dim) + "D");
    }
    if (!(array.flags() & py::array::c_style)) {
        array = py::array_t<T>(array, true);
    }
    return array.request();
}

}  // namespace

PYBIND11_MODULE(_native, m) {
    m.doc() = "Native Helix physics kernel bindings";

    m.def(
        "compute_on_target_score_encoded",
        [](py::array_t<std::uint8_t> guide,
           py::array_t<std::uint8_t> window,
           py::array_t<float> weights,
           float bulge_penalty,
           float gc_low,
           float gc_high,
           float pam_penalty,
           float min_score) {
            auto guide_info = ensure_contiguous(guide, "guide", 1);
            auto window_info = ensure_contiguous(window, "window", 1);
            auto weight_info = ensure_contiguous(weights, "weights", 1);
            PhysicsParams params{
                static_cast<const float*>(weight_info.ptr),
                static_cast<std::size_t>(weight_info.size),
                bulge_penalty,
                gc_low,
                gc_high,
                pam_penalty,
                min_score,
            };
            return compute_on_target_score_encoded(
                static_cast<const std::uint8_t*>(guide_info.ptr),
                static_cast<std::size_t>(guide_info.size),
                static_cast<const std::uint8_t*>(window_info.ptr),
                static_cast<std::size_t>(window_info.size),
                params);
        },
        py::arg("guide"),
        py::arg("window"),
        py::arg("weights"),
        py::arg("bulge_penalty"),
        py::arg("gc_low"),
        py::arg("gc_high"),
        py::arg("pam_penalty"),
        py::arg("min_score"));

    m.def(
        "score_pairs_encoded",
        [](py::array_t<std::uint8_t> guides,
           py::array_t<std::uint8_t> windows,
           py::array_t<float> weights,
           float bulge_penalty,
           float gc_low,
           float gc_high,
           py::object pam_penalties,
           float min_score) {
            auto guide_info = ensure_contiguous(guides, "guides", 2);
            auto window_info = ensure_contiguous(windows, "windows", 2);
            auto weight_info = ensure_contiguous(weights, "weights", 1);
            if (guide_info.shape[1] != window_info.shape[1]) {
                throw py::value_error("guides and windows must share the same sequence length");
            }
            const std::size_t n_guides = static_cast<std::size_t>(guide_info.shape[0]);
            const std::size_t n_windows = static_cast<std::size_t>(window_info.shape[0]);
            const std::size_t seq_len = static_cast<std::size_t>(guide_info.shape[1]);
            const float* pam_ptr = nullptr;
            py::array_t<float> pam_array;
            if (!pam_penalties.is_none()) {
                pam_array = pam_penalties.cast<py::array_t<float>>();
                auto pam_info = ensure_contiguous(pam_array, "pam_penalties", 2);
                if (pam_info.shape[0] != guide_info.shape[0] || pam_info.shape[1] != window_info.shape[0]) {
                    throw py::value_error("pam_penalties shape must match (guides, windows)");
                }
                pam_ptr = static_cast<const float*>(pam_info.ptr);
            }
            PhysicsParams params{
                static_cast<const float*>(weight_info.ptr),
                static_cast<std::size_t>(weight_info.size),
                bulge_penalty,
                gc_low,
                gc_high,
                /*pam_penalty=*/0.0f,
                min_score,
            };
            py::array_t<float> output({guide_info.shape[0], window_info.shape[0]});
            auto out_info = output.request();
            score_pairs_encoded(
                static_cast<const std::uint8_t*>(guide_info.ptr),
                n_guides,
                static_cast<const std::uint8_t*>(window_info.ptr),
                n_windows,
                seq_len,
                params,
                pam_ptr,
                static_cast<float*>(out_info.ptr));
            return output;
        },
        py::arg("guides"),
        py::arg("windows"),
        py::arg("weights"),
        py::arg("bulge_penalty"),
        py::arg("gc_low"),
        py::arg("gc_high"),
        py::arg("pam_penalties") = py::none(),
        py::arg("min_score"));

#ifdef HELIX_ENGINE_HAS_CUDA
    m.def("cuda_available", []() { return helix::engine::cuda::is_available(); });
    m.def(
        "score_pairs_encoded_cuda",
        [](py::array_t<std::uint8_t> guides,
           py::array_t<std::uint8_t> windows,
           py::array_t<float> weights,
           float bulge_penalty,
           float gc_low,
           float gc_high,
           py::array_t<float> pam_penalties,
           float min_score) {
            auto guide_info = ensure_contiguous(guides, "guides", 2);
            auto window_info = ensure_contiguous(windows, "windows", 2);
            auto weight_info = ensure_contiguous(weights, "weights", 1);
            auto pam_info = ensure_contiguous(pam_penalties, "pam_penalties", 2);
            if (guide_info.shape[1] != window_info.shape[1]) {
                throw py::value_error("guides and windows must share the same sequence length");
            }
            if (pam_info.shape[0] != guide_info.shape[0] || pam_info.shape[1] != window_info.shape[0]) {
                throw py::value_error("pam_penalties shape must match (guides, windows)");
            }
            PhysicsParams params{
                static_cast<const float*>(weight_info.ptr),
                static_cast<std::size_t>(weight_info.size),
                bulge_penalty,
                gc_low,
                gc_high,
                0.0f,
                min_score,
            };
            py::array_t<float> output({guide_info.shape[0], window_info.shape[0]});
            auto out_info = output.request();
            helix::engine::cuda::score_pairs_encoded(
                static_cast<const std::uint8_t*>(guide_info.ptr),
                static_cast<std::size_t>(guide_info.shape[0]),
                static_cast<const std::uint8_t*>(window_info.ptr),
                static_cast<std::size_t>(window_info.shape[0]),
                static_cast<std::size_t>(guide_info.shape[1]),
                params,
                static_cast<const float*>(pam_info.ptr),
                static_cast<float*>(out_info.ptr));
            return output;
        },
        py::arg("guides"),
        py::arg("windows"),
        py::arg("weights"),
        py::arg("bulge_penalty"),
        py::arg("gc_low"),
        py::arg("gc_high"),
        py::arg("pam_penalties"),
        py::arg("min_score"));
#else
    m.def("cuda_available", []() { return false; });
#endif
}

}  // namespace helix::engine
