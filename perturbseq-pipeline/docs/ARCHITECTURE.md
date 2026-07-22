# Architecture and General Diagrams

**Author:** Mohamed Kassam  
**Date:** July 19, 2026

## End-to-end workflow

```mermaid
flowchart TD
    A{Starting point} -->|FASTQ| B[FastQC and MultiQC]
    A -->|Cell Ranger H5| E[Input validation]
    B --> C[Cell Ranger count and feature-barcode processing]
    C --> E
    E --> F[Create Seurat object]
    F --> G[Cell and gene QC]
    G --> H[Guide assignment]
    H --> I[Guide confidence and target suppression review]
    I --> J[SCTransform v2]
    J --> K[Cell-cycle scoring / optional regression]
    K --> L[PCA, neighbors, UMAP, clustering]
    L --> M[Metadata-driven model and contrasts]
    M --> N[DGE]
    N --> O[GSEA / ORA]
    O --> P[Tables, RDS, figures, report, provenance]
```

## Metadata-driven contrast model

```mermaid
flowchart LR
    A[Cell metadata] --> D[Design builder]
    B[Reference levels] --> D
    C[Explicit contrast YAML] --> D
    D --> E[Treatment effects]
    D --> F[Perturbation effects]
    D --> G[Conditional effects]
    D --> H[Interaction terms]
```

## Dual-guide decision logic

```mermaid
flowchart TD
    A[Guide counts per cell] --> B{Sufficient guide UMIs?}
    B -->|No| C[NoGuide / LowConfidence]
    B -->|Yes| D{Top-to-second ratio passes?}
    D -->|No| E[Ambiguous / MultiGuide]
    D -->|Yes| F{Guide1 and Guide2 targets}
    F -->|Same gene| G[GeneA/GeneA perturbation]
    F -->|Different genes| H[GeneA/GeneB combination]
    G --> I[Validate target suppression and guide concordance]
    H --> I
```

## Statistical decision framework

```mermaid
flowchart TD
    A{Biological replicates available?} -->|Yes| B[Pseudobulk by sample and condition]
    B --> C[edgeR / limma-voom / DESeq2]
    A -->|No| D[Exploratory cell-level model]
    D --> E[MAST with covariates]
    E --> F[Emphasize effect size and consistency]
    C --> G[Multiple-testing correction]
    F --> G
    G --> H[Ranked enrichment and interpretation]
```

## General example description

This repository separates orchestration, scientific analysis, and experiment-specific configuration. Nextflow provides reproducible execution and provenance; R provides established single-cell and statistical methods; Python performs input validation and design/contrast generation. The main trade-off is flexibility versus strict assumptions: the software remains generic, while the experiment YAML makes controls and desired comparisons explicit and auditable.
