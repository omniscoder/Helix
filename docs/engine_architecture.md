# Helix Engine Architecture

This document captures the now-frozen boundary between Helix's Python shell and
the native execution engines.

## Public API Surface (Python)

- `helix.crispr.simulator` provides dataclasses like `EfficiencyTargetRequest`
  and batch APIs such as `predict_efficiency_for_targets`. This layer is pure
  Python and **must not** grow new responsibilities when chasing performance.
- `helix.crispr.physics` defines the `CrisprPhysics` protocol and helpers such
  as `compute_on_target_score_encoded` / `score_pairs_encoded`. These work
  exclusively with encoded `np.ndarray` buffers and simple numeric parameters.

All callers should treat these as the final, stable boundary. Switching between
`backend="cpu-reference"`, `backend="native-cpu"`, or `backend="gpu"` happens
inside `create_crispr_physics` with no API changes.

### Configuration / CLI Surface

- CLI flag `--engine-backend {cpu-reference,native-cpu,gpu}` (or env
  `HELIX_CRISPR_BACKEND`) sets the preferred backend. Defaults to native-cpu when
  available, otherwise cpu-reference.
- `--engine-allow-fallback` / `HELIX_CRISPR_ALLOW_FALLBACK=1` lets Helix fall
  back to cpu-reference if the requested backend can’t load.
- Studio and other Python callers inherit the same environment variables by
  default; no extra wiring is required beyond the batch APIs above.
- `helix engine info` prints the resolved backend, fallback state, and whether
  the native module is available—handy for bug reports.
- Helix Studio surfaces the active backend/fallback/native status in the control
  panel so users immediately see which engine is powering the visualization.
- All CRISPR artifacts now carry `crispr_scoring_version` (currently 1.0.0) and
  `crispr_engine_backend` metadata. Any deliberate scoring change must bump the
  version, refresh the fixtures, and document the change (see
  `docs/scoring_version.md`).
- Backend behavior matrix (requested backend vs availability):

  | Requested | Native built? | CUDA available? | Fallback allowed? | Result backend | Behavior |
  |-----------|---------------|-----------------|-------------------|----------------|----------|
  | gpu       | yes           | yes             | n/a               | gpu            | CUDA kernels run; metadata reports `gpu`. |
  | gpu       | yes/no        | no              | no                | –              | `RuntimeError` (“Helix CUDA backend is unavailable.”). |
  | gpu       | yes/no        | no              | yes               | cpu-reference  | Warning logged; metadata records `cpu-reference`. |
  | native-cpu| yes           | n/a             | n/a               | native-cpu     | Native extension handles scoring. |
  | native-cpu| no            | n/a             | no                | –              | `RuntimeError` explaining the native build is missing. |
  | native-cpu| no            | n/a             | yes               | cpu-reference  | Warning logged; metadata records `cpu-reference`. |
  | cpu-reference| n/a        | n/a             | n/a               | cpu-reference  | Pure Python reference implementation, always available. |

  `use_gpu=True` in helper APIs simply selects the `gpu` row above.

### GPU Backend (Experimental)

The CUDA path is available when you build `helix_engine._native` with
`HELIX_ENGINE_ENABLE_CUDA=ON`. PyPI wheels remain CPU-only; build locally via:

```
./scripts/build_native_cuda.sh
HELIX_CRISPR_BACKEND=gpu helix engine benchmark --backends gpu --json gpu_benchmark.json
```

If CUDA is unavailable the CLI falls back to `cpu-reference` and records the
fallback in both tables and JSON.

### CRISPR Engine Throughput (Example Dev Machine)

Measured on a CUDA-enabled dev container (RTX A4000, Python 3.12, build `-O3`).
Run `helix engine benchmark` (or `python -m helix.cli engine benchmark`) to
collect the same metrics locally and compare against this table.

| Backend        | G   | N      | L  | MPairs/s |
|----------------|-----|--------|----|---------:|
| cpu-reference  | 1   | 512    | 20 |     0.06 |
| native-cpu     | 1   | 512    | 20 |     0.90 |
| gpu            | 1   | 512    | 20 |    <0.01 |
| cpu-reference  | 96  | 4096   | 20 |     0.07 |
| native-cpu     | 96  | 4096   | 20 |     7.53 |
| gpu            | 96  | 4096   | 20 |     3.31 |
| cpu-reference  | 256 | 16384  | 20 |     0.09 |
| native-cpu     | 256 | 16384  | 20 |    11.05 |
| gpu            | 256 | 16384  | 20 |    12.17 |

Tiny GPU workloads (1×512×20) are dominated by kernel launch overhead and settle around 0.002 MPairs/s even though they quantize to 0.00 at two decimal places. Larger batches amortize the cost and reach parity/speedups vs native.

Prime editing respects the same knobs: by default it mirrors the CRISPR backend choice, but you can override it via `HELIX_PRIME_BACKEND` and `HELIX_PRIME_ALLOW_FALLBACK` if needed. Prime physics scoring terminology is documented in `docs/prime_physics.md`.

### Prime Engine Throughput (Example Dev Machine)

| Backend       | Targets | Genome nt | Spacer | Predictions/sec |
|---------------|---------|-----------|--------|----------------:|
| cpu-reference | 32      | 2,000     | 20     |          5634.48 |
| cpu-reference | 64      | 4,000     | 20     |          3037.36 |

### Benchmark JSON Schema (v1)

`helix engine benchmark --json` emits a structured payload so CI/doc tooling can
track regressions. Schema v1 fields:

- `helix_version` – package version string.
- `scoring_versions.crispr` / `.prime` – scoring-version identifiers.
- `env.platform`, `env.python_version`, `env.cuda_available`, `env.gpu_name`,
  `env.native_backend_available`, `env.backends_built`.
- `seed` – RNG seed used for the synthetic data.
- `config.backends`, `config.crispr_shapes`, `config.prime_workloads`.
- `benchmarks.crispr[]` entries contain `backend_requested`, `backend_used`,
  `shape` (`GxNxL`), raw `g/n/l`, `mpairs_per_s`, and `elapsed_seconds`.
- `benchmarks.prime[]` entries contain `backend_requested`, `backend_used`,
  `workload` (`targetsxgenome_lenxspacer`), `predictions`/`predictions_per_s`,
  and `elapsed_seconds`.

Example (trimmed) JSON:

```json
{
  "helix_version": "0.4.0",
  "scoring_versions": {"crispr": "1.0.0", "prime": "1.0.0"},
  "env": {
    "platform": "Linux-6.6.87-x86_64",
    "python_version": "3.12.3",
    "cuda_available": false,
    "gpu_name": null,
    "native_backend_available": false,
    "backends_built": ["cpu-reference"]
  },
  "seed": 1,
  "config": {
    "backends": ["cpu-reference"],
    "crispr_shapes": ["1x512x20"],
    "prime_workloads": ["32x2000x20"]
  },
  "benchmarks": {
    "crispr": [
      {
        "backend_requested": "cpu-reference",
        "backend_used": "cpu-reference",
        "shape": "1x512x20",
        "mpairs_per_s": 0.06,
        "elapsed_seconds": 0.53
      }
    ],
    "prime": [
      {
        "workload": "32x2000x20",
        "predictions": 32,
        "predictions_per_s": 5634.48,
        "elapsed_seconds": 0.006
      }
    ]
  }
}
```
Future versions of the schema will only add fields; existing keys remain stable.

## Remote Engine Contract

Remote runners (OGN) reuse the same JSON contracts as the CLI:

- Performance endpoints must emit the benchmark schema above.
- Prime scoring endpoints must include the `physics_score` block defined in
  `docs/prime_physics.md`.

See `docs/ogn_integration.md` for the exact fields expected from remote
benchmarks and scoring services.

## Kernel Contract

Every backend (Python, native C++, CUDA) implements the same kernel semantics:

- Guides are encoded as uint8 arrays with shape `(G, L)` (row-major).
- Windows are encoded as uint8 arrays with shape `(N, L)`.
- Output scores are float32 arrays with shape `(G, N)`.
- No broadcasting is permitted; mismatched shapes raise errors.
- Encodings follow the `helix.crispr.physics._encode_sequence_to_uint8`
  convention: A/C/G/T → 0/1/2/3 with all other bases mapping to 4.
- The shared helpers live in `helix.engine.encoding`, so CRISPR, Prime, and any
  other engines use the exact same mapping.

This is the exact contract exposed via the C++ header in
`src/helix_engine/include/helix/physics.hpp` and mirrored by CUDA kernels in the
future.

## Native Engine Layout

```
src/helix_engine/
├── include/helix/physics.hpp   # PhysicsParams + C++ kernel prototypes
├── src/physics.cpp             # Reference CPU implementations
├── src/binding.cpp             # pybind11 bindings (exposes _native module)
└── native.py                   # Python shim that calls into the extension
```

`helix_engine.native` exposes the same scoring functions to Python and only
falls back to pure Python helpers when the extension is missing. When users
explicitly request `backend="native-cpu"`, Helix requires the compiled module
and raises a clear error otherwise.

## Testing Expectations

- `tests/test_crispr_simulator.py` keeps exercising the high-level simulator and
  ensures encoded helpers behave consistently.
- `tests/test_crispr_native_engine.py` (skipped unless `_native` is importable)
  asserts scalar/batch parity between the C++ backend and the Python reference.
- Golden fixtures live in `tests/data/engine_crispr`. The test suite compares
  every backend (`cpu-reference`, `native-cpu`, later `gpu`) against those saved
  scores (`tests/test_crispr_engine_fixtures.py`). Any change to scoring must
  update the fixture data explicitly.
- Prime fixtures live alongside them in `tests/data/engine_prime`, enforced by
  `tests/test_prime_engine.py`.
- CI environments without the native build remain fully green; local developer
  environments with the extension compiled gain stronger guarantees.

## Performance Harness

- `benchmarks/engine_crispr_perf.py` generates random encoded arrays and
  measures throughput for `score_pairs_encoded` across user-specified `(G, N, L)`
  regimes and backends. Run it locally to track wins/regressions before landing
  CUDA or other optimizations. Example (reference backend on dev container):

  ```
  backend=cpu-reference shape=(1,512,20)   -> 0.05 MPairs/s
  backend=cpu-reference shape=(32,4096,20) -> 0.06 MPairs/s
  ```

- For comparing multiple backends/configs, call the helpers in
  `helix.crispr.physics` (e.g. `score_pairs_encoded_multi`) at a higher layer
  rather than overloading the core kernel API. The kernel itself always assumes
  “one backend, many guides/windows.”
 - Prime editing follows the same pattern. `benchmarks/engine_prime_perf.py`
  offers a starter perf harness, and `predict_prime_outcomes_for_targets`
  provides the batch API until the physics engine matures further. Prime JSON
  outputs report `prime_scoring_version` (currently 1.0.0) and
  `prime_engine_backend` fields so consumers know which engine generated them.
  Example baseline (`targets=64`, `genome_len=2000`, spacer length 20) on this
  dev container: `cpu-reference` ≈ 5.4K predictions/sec.

With this structure in place, adding CUDA support means implementing the same
`score_pairs_encoded` kernel on the device and hooking it up to the existing
protocol—no changes to simulator or benchmark code are required. A native CUDA
backend is available when the engine is built with
`-DHELIX_ENGINE_ENABLE_CUDA=ON`; see `docs/engine_cuda_plan.md` for build notes
and future milestones.
