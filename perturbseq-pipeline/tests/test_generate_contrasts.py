"""Unit tests for metadata-driven contrast generation.

Author: Mohamed Kassam
Date: 2026-07-19
"""
import importlib.util
from pathlib import Path
import pandas as pd

MODULE_PATH = Path(__file__).parents[1] / "bin" / "generate_contrasts.py"
spec = importlib.util.spec_from_file_location("generate_contrasts", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_custom_reference_names_are_supported():
    metadata = pd.DataFrame({
        "condition": ["Untreated", "CompoundX", "Untreated", "CompoundX"],
        "target": ["SafeGuide", "SafeGuide", "GeneA", "GeneA"],
    })
    config = {
        "metadata": {"treatment_column": "condition", "perturbation_column": "target"},
        "references": {"treatment": "Untreated", "perturbation": "SafeGuide"},
        "contrasts": {"auto_generate": {"treatment_vs_reference": True, "perturbation_vs_reference": True, "interaction": True}, "explicit": []},
    }
    result = module.generate(metadata, config)
    assert "CompoundX_vs_Untreated" in set(result.contrast_id)
    assert "GeneA_vs_SafeGuide" in set(result.contrast_id)
    assert "CompoundX_x_GeneA" in set(result.contrast_id)
    assert not result.astype(str).apply(lambda col: col.str.contains("DMSO", regex=False)).any().any()
