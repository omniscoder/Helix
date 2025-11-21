# GPU Profiling Cheatsheet

Use these commands when you need deep visibility into the CUDA kernels. They
assume you already built `_native` with CUDA (`./scripts/build_native_cuda.sh`).

## Nsight Systems (timeline)

```
nsys profile -o helix_gpu_timeline \
  env HELIX_CRISPR_BACKEND=gpu PYTHONPATH=src:. \
  python -m helix.cli engine benchmark --backends gpu --crispr-shapes 256x16384x20
```

This captures kernel launches, CPU/GPU overlap, and overall timeline. Open the
`.qdrep` file in Nsight Systems to inspect occupancy/stalls.

## Nsight Compute (kernel metrics)

```
ncu -o helix_gpu_kernel \
  env HELIX_CRISPR_BACKEND=gpu PYTHONPATH=src:. \
  python -m helix.cli engine benchmark --backends gpu --crispr-shapes 96x4096x20
```

Look at:

- Achieved occupancy
- DRAM throughput vs the GPU’s peak
- Warp execution efficiency (branch divergence)
- SM utilization

These metrics highlight whether optimization should focus on memory, math, or
branch divergence. Use smaller shapes for quick iteration and large shapes for
steady-state analysis.
