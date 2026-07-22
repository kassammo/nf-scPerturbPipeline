#!/usr/bin/env python3
# Author: Mohamed Kassam
# Date: 2026-07-19
"""Metadata-driven Perturb-seq matrix workflow.

This lightweight integration runner is intended for smoke/integration testing and
small datasets. Production statistical inference should use the R/Seurat-MAST or
pseudobulk implementation included in the repository when biological replicates
are available.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from sklearn.decomposition import PCA


def bh_adjust(pvalues: np.ndarray | pd.Series) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def read_config(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise ValueError("Experiment config must contain a YAML mapping.")
    return cfg


def validate_inputs(counts: pd.DataFrame, guides: pd.DataFrame, meta: pd.DataFrame, cfg: dict[str, Any]) -> None:
    barcode_col = cfg["metadata"]["barcode_column"]
    missing = {barcode_col} - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")
    if counts.columns.duplicated().any() or guides.columns.duplicated().any():
        raise ValueError("Cell barcodes must be unique in count matrices.")
    if not set(counts.columns).issubset(set(meta[barcode_col])):
        raise ValueError("Every expression-matrix barcode must exist in metadata.")
    if set(counts.columns) != set(guides.columns):
        raise ValueError("Gene and guide matrices must contain the same cell barcodes.")
    for factor in cfg["experiment"].get("factors", []):
        column = factor["column"]
        if column not in meta.columns:
            raise ValueError(f"Configured factor column '{column}' is absent from metadata.")
        reference = factor.get("reference")
        if reference is not None and reference not in set(meta[column].astype(str)):
            raise ValueError(f"Reference level '{reference}' is absent from metadata column '{column}'.")


def assign_guides(guides: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    guide_cfg = cfg.get("guides", {})
    min_umi = int(guide_cfg.get("min_umi", 3))
    min_ratio = float(guide_cfg.get("min_ratio_to_second", 2.0))
    no_guide = str(guide_cfg.get("no_guide_label", "NoGuide"))
    target_regex = str(guide_cfg.get("target_extraction_regex", r"^(.+?)(?:[_-]g?\d+)?$"))
    non_targeting_regex = str(guide_cfg.get("non_targeting_regex", r"(?i)NTC|non[-_ ]?target|safe"))
    non_targeting_label = str(guide_cfg.get("non_targeting_label", "NTC"))

    arr = guides.to_numpy()
    top_two = np.argsort(-arr, axis=0)[:2]
    top = arr[top_two[0], np.arange(arr.shape[1])]
    second = arr[top_two[1], np.arange(arr.shape[1])] if arr.shape[0] > 1 else np.zeros(arr.shape[1])
    top_guides = guides.index.to_numpy()[top_two[0]]
    confident = (top >= min_umi) & (top / np.maximum(second, 1) >= min_ratio)

    targets: list[str] = []
    for guide, ok in zip(top_guides, confident):
        if not ok:
            targets.append(no_guide)
        elif re.search(non_targeting_regex, str(guide)):
            targets.append(non_targeting_label)
        else:
            match = re.search(target_regex, str(guide))
            targets.append(match.group(1) if match else str(guide))

    return pd.DataFrame(
        {
            "assigned_guide": np.where(confident, top_guides, no_guide),
            "assigned_perturbation": targets,
            "guide_top_umi": top,
            "guide_second_umi": second,
            "guide_ratio": top / np.maximum(second, 1),
            "guide_confident": confident,
        },
        index=guides.columns,
    )


def build_group_labels(meta: pd.DataFrame, cfg: dict[str, Any]) -> pd.Series:
    columns = cfg["experiment"].get("group_columns")
    if not columns:
        columns = [factor["column"] for factor in cfg["experiment"].get("factors", [])]
    separator = str(cfg["experiment"].get("group_separator", "__"))
    for column in columns:
        if column not in meta.columns:
            raise ValueError(f"Group column '{column}' is absent from metadata.")
    return meta[columns].astype(str).agg(separator.join, axis=1)


def generate_contrasts(meta: pd.DataFrame, cfg: dict[str, Any]) -> list[dict[str, str]]:
    contrast_cfg = cfg.get("contrasts", {})
    explicit = contrast_cfg.get("explicit", []) or []
    contrasts: list[dict[str, str]] = []
    for item in explicit:
        contrasts.append({"name": str(item["name"]), "group1": str(item["group1"]), "group0": str(item["group0"])})

    if contrast_cfg.get("auto_all_vs_reference", False):
        reference = contrast_cfg.get("reference_group")
        if not reference:
            reference_values = []
            factor_by_column = {f["column"]: f for f in cfg["experiment"].get("factors", [])}
            for column in cfg["experiment"].get("group_columns", factor_by_column):
                factor = factor_by_column.get(column, {})
                if "reference" not in factor:
                    raise ValueError(f"No reference level configured for group factor '{column}'.")
                reference_values.append(str(factor["reference"]))
            reference = str(cfg["experiment"].get("group_separator", "__")).join(reference_values)
        groups = sorted(set(meta["analysis_group"].astype(str)))
        for group in groups:
            if group != reference:
                contrasts.append({"name": f"{group}_vs_{reference}", "group1": group, "group0": str(reference)})

    seen: set[str] = set()
    unique: list[dict[str, str]] = []
    for contrast in contrasts:
        if contrast["name"] not in seen:
            unique.append(contrast)
            seen.add(contrast["name"])
    return unique


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def run() -> None:
    parser = argparse.ArgumentParser(description="Run a metadata-driven Perturb-seq matrix workflow.")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    out = Path(args.outdir)
    (out / "dge").mkdir(parents=True, exist_ok=True)
    (out / "pathway").mkdir(parents=True, exist_ok=True)
    cfg = read_config(Path(args.config))

    files_cfg = cfg.get("files", {})
    counts = pd.read_csv(input_dir / files_cfg.get("gene_counts", "gene_counts.tsv"), sep="\t", index_col=0)
    guides = pd.read_csv(input_dir / files_cfg.get("guide_counts", "guide_counts.tsv"), sep="\t", index_col=0)
    meta = pd.read_csv(input_dir / files_cfg.get("cell_metadata", "cell_metadata.tsv"), sep="\t")
    barcode_col = cfg["metadata"]["barcode_column"]
    validate_inputs(counts, guides, meta, cfg)
    meta = meta.set_index(barcode_col).loc[counts.columns].copy()

    guide_assignment = assign_guides(guides, cfg)
    for column in guide_assignment.columns:
        meta[column] = guide_assignment[column]

    perturbation_column = cfg["metadata"].get("perturbation_column")
    if perturbation_column:
        if perturbation_column not in meta.columns:
            raise ValueError(f"Configured perturbation column '{perturbation_column}' is absent from metadata.")
    else:
        perturbation_column = "assigned_perturbation"
    meta["analysis_perturbation"] = meta[perturbation_column].astype(str)

    total = counts.sum(axis=0)
    features = (counts > 0).sum(axis=0)
    mito_regex = str(cfg.get("qc", {}).get("mitochondrial_regex", r"^MT-"))
    mito_mask = counts.index.to_series().str.contains(mito_regex, regex=True)
    percent_mito = counts.loc[mito_mask].sum(axis=0) / np.maximum(total, 1) * 100
    meta["nCount_RNA"] = total
    meta["nFeature_RNA"] = features
    meta["percent_mito"] = percent_mito

    qc = cfg.get("qc", {})
    keep = (
        (features >= int(qc.get("min_features", 50)))
        & (features <= int(qc.get("max_features", counts.shape[0])))
        & (percent_mito <= float(qc.get("max_percent_mito", 25)))
    )
    if bool(qc.get("require_confident_guide", True)):
        keep &= meta["guide_confident"].astype(bool)
    meta["qc_pass"] = keep

    filtered_counts = counts.loc[:, keep]
    filtered_meta = meta.loc[keep].copy()
    filtered_meta["analysis_group"] = build_group_labels(filtered_meta, cfg)

    normalization_cfg = cfg.get("normalization", {})
    scale_factor = float(normalization_cfg.get("scale_factor", 10000))
    normalized = np.log1p(filtered_counts.div(filtered_counts.sum(axis=0), axis=1) * scale_factor)
    normalized.to_csv(out / "normalized_matrix.tsv", sep="\t")
    filtered_meta.to_csv(out / "cell_metadata.tsv", sep="\t")

    n_components = min(int(cfg.get("reduction", {}).get("npcs", 10)), normalized.shape[1] - 1, normalized.shape[0])
    pca = PCA(n_components=max(2, n_components), random_state=int(cfg.get("seed", 42))).fit_transform(normalized.T.to_numpy())
    pca_df = pd.DataFrame(pca, index=normalized.columns, columns=[f"PC{i + 1}" for i in range(pca.shape[1])])
    pca_df.to_csv(out / "pca.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(8, 6))
    for label, indices in filtered_meta.groupby("analysis_group").groups.items():
        ax.scatter(pca_df.loc[indices, "PC1"], pca_df.loc[indices, "PC2"], s=16, label=label)
    ax.set(xlabel="PC1", ylabel="PC2", title="Perturb-seq PCA by configured analysis group")
    ax.legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "pca_groups.png", dpi=150)
    plt.close(fig)

    contrasts = generate_contrasts(filtered_meta, cfg)
    minimum_cells = int(cfg.get("dge", {}).get("minimum_cells_per_group", 3))
    summaries: list[dict[str, Any]] = []
    for contrast in contrasts:
        group1, group0 = contrast["group1"], contrast["group0"]
        cells1 = filtered_meta.index[filtered_meta["analysis_group"] == group1]
        cells0 = filtered_meta.index[filtered_meta["analysis_group"] == group0]
        if len(cells1) < minimum_cells or len(cells0) < minimum_cells:
            summaries.append({**contrast, "status": "skipped", "reason": "insufficient_cells", "n_group1": len(cells1), "n_group0": len(cells0)})
            continue
        mean1 = normalized[cells1].mean(axis=1)
        mean0 = normalized[cells0].mean(axis=1)
        t_stat, p_values = stats.ttest_ind(
            normalized[cells1].to_numpy(), normalized[cells0].to_numpy(), axis=1, equal_var=False, nan_policy="omit"
        )
        result = pd.DataFrame(
            {
                "gene": normalized.index,
                "avg_log2FC": (mean1 - mean0).to_numpy() / np.log(2),
                "statistic": t_stat,
                "p_value": p_values,
                "p_adj": bh_adjust(p_values),
                "mean_group1": mean1.to_numpy(),
                "mean_group0": mean0.to_numpy(),
                "group1": group1,
                "group0": group0,
                "contrast": contrast["name"],
            }
        ).sort_values("p_adj")
        filename = safe_filename(contrast["name"]) + ".tsv"
        result.to_csv(out / "dge" / filename, sep="\t", index=False)
        summaries.append(
            {
                **contrast,
                "status": "completed",
                "n_group1": len(cells1),
                "n_group0": len(cells0),
                "n_significant_fdr_0.05": int((result["p_adj"] < 0.05).sum()),
                "file": filename,
            }
        )

    gene_set_file = input_dir / files_cfg.get("gene_sets", "gene_sets.tsv")
    pathway_rows: list[dict[str, Any]] = []
    if cfg.get("pathway", {}).get("enabled", True) and gene_set_file.exists():
        gene_sets = pd.read_csv(gene_set_file, sep="\t")
        universe = set(normalized.index)
        fdr_cutoff = float(cfg.get("pathway", {}).get("dge_fdr", 0.1))
        logfc_cutoff = float(cfg.get("pathway", {}).get("dge_abs_log2fc", 0.5))
        for summary in summaries:
            if summary.get("status") != "completed":
                continue
            de = pd.read_csv(out / "dge" / summary["file"], sep="\t")
            significant = set(de.loc[(de["p_adj"] < fdr_cutoff) & (de["avg_log2FC"].abs() > logfc_cutoff), "gene"])
            for set_name, group in gene_sets.groupby("set"):
                genes = set(group["gene"]) & universe
                overlap = significant & genes
                p_value = stats.hypergeom.sf(len(overlap) - 1, len(universe), len(genes), len(significant)) if significant else 1.0
                pathway_rows.append(
                    {
                        "contrast": summary["name"],
                        "pathway": set_name,
                        "overlap": len(overlap),
                        "set_size": len(genes),
                        "p_value": p_value,
                        "genes": ",".join(sorted(overlap)),
                    }
                )
    pathway_df = pd.DataFrame(pathway_rows)
    if len(pathway_df):
        pathway_df["p_adj"] = pathway_df.groupby("contrast")["p_value"].transform(lambda x: bh_adjust(x.to_numpy()))
    pathway_df.to_csv(out / "pathway" / "ora.tsv", sep="\t", index=False)

    expected_col = cfg["metadata"].get("expected_perturbation_column")
    accuracy = None
    if expected_col and expected_col in filtered_meta.columns:
        accuracy = float((filtered_meta["analysis_perturbation"] == filtered_meta[expected_col].astype(str)).mean())

    summary = {
        "input_cells": int(counts.shape[1]),
        "retained_cells": int(filtered_counts.shape[1]),
        "genes": int(filtered_counts.shape[0]),
        "analysis_groups": filtered_meta["analysis_group"].value_counts().to_dict(),
        "guide_assignment_accuracy_vs_expected": accuracy,
        "contrasts": summaries,
        "config": str(Path(args.config)),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    html = (
        "<html><body><h1>Metadata-driven Perturb-seq Report</h1>"
        f"<pre>{json.dumps(summary, indent=2)}</pre>"
        "<img src='pca_groups.png' width='900'>"
        "<p>Outputs: normalized_matrix.tsv, cell_metadata.tsv, pca.tsv, dge/, pathway/.</p>"
        "</body></html>"
    )
    (out / "report.html").write_text(html)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    run()
