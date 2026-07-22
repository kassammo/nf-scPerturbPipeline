# Example Presentation Guide

**Author:** Mohamed Kassam  
**Date:** July 19, 2026

## Two-minute overview

The repository implements a modular Perturb-seq workflow from Cell Ranger outputs to guide QC, normalization, dimensionality reduction, differential expression, interaction analysis, pathway enrichment, and reproducible reporting. The main engineering choice was to keep all biological labels and reference conditions outside the source code. A user defines control levels and contrasts in YAML, allowing the same pipeline to support vehicle controls, untreated controls, stimulation experiments, single perturbations, and combinatorial designs.

## Key design choices to explain

- Nextflow DSL2 for portability, provenance, retries, and resource control.
- Seurat/SCTransform for established single-cell preprocessing.
- Explicit guide-confidence rules rather than accepting every detected barcode.
- Replicate-aware pseudobulk inference when possible; MAST for exploratory cell-level analysis.
- Metadata-driven contrasts to prevent DMSO- or gene-specific hard-coding.
- Cell Ranger kept external because of licensing and reference-size constraints.

## Trade-offs

Cell-cycle regression can remove confounding but may also remove real perturbation biology. Strict guide thresholds improve specificity but reduce cell yield. Cell-level methods provide sensitivity but risk pseudoreplication. Automatic contrasts improve usability, while explicit contrasts remain necessary for scientifically precise questions.

## Honest validation statement

The repository includes software unit tests and CI checks. It should not be presented as biologically validated until run against a representative public or client-approved Perturb-seq dataset with predefined acceptance criteria.
