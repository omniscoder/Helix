# Contributing

1. Fork → `pip install -e ".[dev]"`.
2. Run `pytest -q` plus `mkdocs serve` to preview docs.
3. Open a PR with:
   - Tests for new features (especially schema/viz contracts).
   - Docs updates (README + MkDocs page).
   - `CHANGELOG.md` entry.

Schema or viz changes should include manifest diffs and sample artifacts (`helix demo viz`) so reviewers can verify provenance.
