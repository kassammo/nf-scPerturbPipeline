#!/usr/bin/env python3
"""Generate metadata-driven contrasts from an experiment configuration.

Author: Mohamed Kassam
Date: 2026-07-19
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import yaml


def _filter_expression(assignments: dict[str, object]) -> str:
    return " and ".join(f'`{column}` == {value!r}' for column, value in assignments.items())


def generate(metadata: pd.DataFrame, config: dict) -> pd.DataFrame:
    meta_cfg = config.get("metadata", {})
    refs = config.get("references", {})
    contrast_cfg = config.get("contrasts", {})
    auto = contrast_cfg.get("auto_generate", {})
    treatment_col = meta_cfg.get("treatment_column", "treatment")
    perturbation_col = meta_cfg.get("perturbation_column", "perturbation")
    treatment_ref = refs.get("treatment")
    perturbation_ref = refs.get("perturbation")

    required = {treatment_col, perturbation_col}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing configured columns: {sorted(missing)}")

    rows: list[dict[str, str]] = []
    treatments = sorted(metadata[treatment_col].dropna().astype(str).unique())
    perturbations = sorted(metadata[perturbation_col].dropna().astype(str).unique())

    if auto.get("treatment_vs_reference", True) and treatment_ref is not None:
        for treatment in treatments:
            if treatment == str(treatment_ref):
                continue
            numerator = {treatment_col: treatment, perturbation_col: perturbation_ref}
            denominator = {treatment_col: treatment_ref, perturbation_col: perturbation_ref}
            rows.append({"contrast_id": f"{treatment}_vs_{treatment_ref}", "group1_filter": _filter_expression(numerator), "group0_filter": _filter_expression(denominator), "contrast_type": "treatment"})

    if auto.get("perturbation_vs_reference", True) and perturbation_ref is not None:
        for perturbation in perturbations:
            if perturbation == str(perturbation_ref):
                continue
            numerator = {treatment_col: treatment_ref, perturbation_col: perturbation}
            denominator = {treatment_col: treatment_ref, perturbation_col: perturbation_ref}
            rows.append({"contrast_id": f"{perturbation}_vs_{perturbation_ref}", "group1_filter": _filter_expression(numerator), "group0_filter": _filter_expression(denominator), "contrast_type": "perturbation"})

    if auto.get("interaction", True) and treatment_ref is not None and perturbation_ref is not None:
        for treatment in treatments:
            if treatment == str(treatment_ref):
                continue
            for perturbation in perturbations:
                if perturbation == str(perturbation_ref):
                    continue
                rows.append({"contrast_id": f"{treatment}_x_{perturbation}", "group1_filter": _filter_expression({treatment_col:treatment, perturbation_col:perturbation}), "group0_filter": _filter_expression({treatment_col:treatment, perturbation_col:perturbation_ref}), "contrast_type": "conditional_perturbation"})

    for item in contrast_cfg.get("explicit", []):
        rows.append({"contrast_id": item["id"], "group1_filter": _filter_expression(item["numerator"]), "group0_filter": _filter_expression(item["denominator"]), "contrast_type": "explicit"})

    return pd.DataFrame(rows).drop_duplicates(subset=["contrast_id"]) if rows else pd.DataFrame(columns=["contrast_id","group1_filter","group0_filter","contrast_type"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--experiment", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    with args.experiment.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    result = generate(pd.read_csv(args.metadata), config)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)
    print(f"Generated {len(result)} contrasts: {args.out}")


if __name__ == "__main__":
    main()
