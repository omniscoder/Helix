# CUDA Backend Plan

The CRISPR engine surface is frozen as `(G, L) × (N, L) → (G, N)`. The first
CUDA milestone should:

1. Implement a `score_pairs_encoded` kernel in CUDA (`guides[G,L]`,
   `windows[N,L]`, `scores[G,N]`).
2. Mirror the existing C++ CPU semantics exactly (same encoding, penalties,
   scoring version).
3. Expose it via `backend="gpu"` in `create_crispr_physics`.
4. Validate against the golden fixtures (`tests/data/engine_crispr/`).
5. Benchmark with `benchmarks/engine_crispr_perf.py --backend gpu ...`.

Keep GPU-specific padding or alignment hidden inside the backend; Python and the
simulator should never see it.

## Building with CUDA

`src/helix_engine/CMakeLists.txt` exposes `HELIX_ENGINE_ENABLE_CUDA=ON`. When
set, CMake enables the CUDA language and compiles `src/physics_cuda.cu` (currently
just a stub). Example:

```bash
cmake -S src/helix_engine -B build/helix_engine -DHELIX_ENGINE_ENABLE_CUDA=ON
```

Future patches will replace the stub with a real kernel and plumb the Python
backend selection through `backend="gpu"`.
