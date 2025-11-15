from __future__ import annotations

import textwrap

from pathlib import Path

from helix.veribiota import dag_payload_to_lean, dag_payloads_to_lean, module_name_from_path


def _sample_payload() -> dict:
    return {
        "nodes": {
            "root": {
                "log_prob": 0.0,
                "metadata": {"stage": "root"},
                "parent_ids": [],
                "seq_hashes": {"chr1": "abc"},
                "diffs": [],
                "sequences": {"chr1": "ACGT"},
            },
            "child": {
                "log_prob": -0.5,
                "metadata": {"stage": "cut"},
                "parent_ids": ["root"],
                "seq_hashes": {"chr1": "def"},
                "diffs": [
                    {
                        "chrom": "chr1",
                        "start": 1,
                        "end": 2,
                        "replacement": "C",
                        "metadata": {"label": "clean_cut"},
                    }
                ],
                "sequences": {"chr1": "AGCT"},
            },
        },
        "edges": [
            {
                "source": "root",
                "target": "child",
                "rule": "crispr.clean_cut",
                "event": {
                    "chrom": "chr1",
                    "start": 1,
                    "end": 1,
                    "replacement": "",
                    "metadata": {"label": "clean_cut"},
                },
            }
        ],
        "root_id": "root",
    }


def test_dag_payload_to_lean_nested_namespace() -> None:
    payload = _sample_payload()
    lean_text = dag_payload_to_lean(
        payload,
        dag_name="microDag",
        module_name="VeriBiota.Bridge",
        import_module="VeriBiota",
        include_eval=False,
        include_theorem=False,
    )
    assert "namespace VeriBiota" in lean_text
    assert "namespace Bridge" in lean_text
    assert "def microDag : EditDAG :=" in lean_text
    assert "rootId := \"root\"" in lean_text
    assert "logProb := -0.5" in lean_text
    assert "#eval" not in lean_text
    expected_tail = textwrap.dedent(
        """
        end Bridge

        end VeriBiota
        """
    ).strip()
    assert expected_tail in lean_text


def test_dag_payloads_to_lean_multi() -> None:
    payload = _sample_payload()
    lean_text = dag_payloads_to_lean(
        [payload, payload],
        ["example_one", "example_two"],
        module_name="Helix.CrisprExamples",
        list_name="exampleDags",
        theorem_name="examples_all_checked",
        include_eval=False,
        include_theorem=True,
    )
    assert "namespace Helix" in lean_text
    assert "namespace CrisprExamples" in lean_text
    assert "def example_one : EditDAG :=" in lean_text
    assert "def example_two : EditDAG :=" in lean_text
    assert "def exampleDags : List EditDAG :=" in lean_text
    assert "theorem examples_all_checked" in lean_text
    assert "∀ dag ∈ exampleDags" in lean_text


def test_module_name_from_path() -> None:
    module_path = Path("Biosim/VeriBiota/Helix/CrisprMicro.lean")
    assert module_name_from_path(module_path) == "Biosim.VeriBiota.Helix.CrisprMicro"
