# Helix 0.4.0 Release Checklist

Follow this checklist every time you cut v0.4.0 (or use it as a template for future releases).

1. Build & Upload

```bash
python -m build
# Inspect dist/
twine upload dist/*
```

2. Tag & Push

```bash
git tag -a v0.4.0 -m "Helix 0.4.0 – engine rewrite + CUDA/Prime/benchmark"
git push origin v0.4.0
```

3. Fresh Environment Verification

```bash
python -m venv /tmp/helix-040-test
source /tmp/helix-040-test/bin/activate
pip install veri-helix
helix engine info
helix engine benchmark --json /tmp/helix-040-benchmark.json
```

Inspect `/tmp/helix-040-benchmark.json` and compare MPairs/s + preds/sec with
`docs/engine_architecture.md` tables. Keep the virtual environment (or archive
the JSON) for release notes.

4. Convenience Script

Use `scripts/release_0_4_0.sh` to run steps 1–2 automatically. Always perform
the clean-environment verification afterward.
