# Pipeline Understanding and Design Rationale

**Author:** Mohamed Kassam  
**Date:** July 19, 2026

## Scientific objective

The pipeline is designed to identify transcriptional programs attributable to treatment, CRISPR perturbation, and their combined effect in Perturb-seq experiments. It separates technical processing from biological configuration so the same implementation can be reused across drugs, stimuli, genes, cell lines, doses, and timepoints.

## Core design principle

Treatment and perturbation names are data, not code. Reference levels and comparisons are defined in an experiment YAML supplied by the runner. The analysis implementation therefore does not assume DMSO, NTC, rapamycin, MTOR, or any other specific label.

## Guide assignment

Guide calls are based on guide abundance and confidence thresholds. Cells can be labelled as no-guide, non-targeting, confidently perturbed, multi-guide, or low-confidence. Where possible, guide assignments should be supported by target-expression suppression and concordance across independent guides.

## Normalization and cell cycle

SCTransform v2 is the default normalization strategy. Cell-cycle scores may be included as covariates when cycle-driven variation obscures the biological question. Regression is not automatic scientific truth: if the perturbation genuinely changes proliferation, aggressive cell-cycle regression could remove meaningful biology. The decision should be documented per experiment.

## Differential expression

The workflow supports cell-level MAST analysis for exploratory work and should use replicate-aware pseudobulk methods for confirmatory inference when biological replicates exist. Drug-by-perturbation questions are represented through conditional comparisons and interaction terms.

## Double perturbations

For guide pairs targeting the same gene, the pair can be collapsed to a gene-level perturbation after confidence and concordance checks. For future GeneA/GeneB designs, the pair should be preserved as an ordered or canonical combination and modeled separately from either single perturbation.

## Pathway analysis

Rank-based enrichment is preferred because it uses the full ordered gene list and avoids a single significance cutoff. ORA can be included as a complementary summary. Gene-set collection versions and organism mappings should be recorded.

## Assumptions

- Cell barcodes can be matched consistently across expression, guide, and metadata tables.
- Controls are represented in the experiment metadata.
- Guide-capture counts are available or perturbation labels have already been assigned.
- Biological interpretation considers experimental replication and guide performance.

## Production-readiness considerations

Before deployment, validate the pipeline on a representative public or client-approved dataset, freeze software and gene-set versions, document acceptance criteria, add dataset-level regression tests, and review computational resources on the target HPC or cloud platform.
