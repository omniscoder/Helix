#!/usr/bin/env bash
set -euo pipefail

echo "==> Building sdist/wheel"
python -m build

echo "==> Uploading to PyPI"
twine upload dist/*

if ! git rev-parse -q --verify refs/tags/v0.4.0 >/dev/null; then
  echo "==> Tagging release"
  git tag -a v0.4.0 -m "Helix 0.4.0 – engine rewrite + CUDA/Prime/benchmark"
fi

echo "==> Pushing tag"
git push origin v0.4.0
