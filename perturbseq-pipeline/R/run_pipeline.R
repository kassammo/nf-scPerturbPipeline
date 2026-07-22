#!/usr/bin/env Rscript
# Author: Mohamed Kassam
# Date: 2026-07-19

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option("--samplesheet", type = "character"),
  make_option("--config", type = "character"),
  make_option("--outdir", type = "character")
)
o <- parse_args(OptionParser(option_list = option_list))
if (is.null(o$samplesheet) || is.null(o$config) || is.null(o$outdir)) {
  stop("--samplesheet, --config, and --outdir are required")
}

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("R/run_pipeline.R")
project_dir <- dirname(dirname(script_path))
source(file.path(project_dir, "R", "utils.R"))
source(file.path(project_dir, "R", "steps.R"))

cfg <- read_yaml(o$config)
set.seed(cfg$seed %||% 42)
ss <- fread(o$samplesheet)
safe_dir(o$outdir)

objs <- lapply(seq_len(nrow(ss)), function(i) import_sample(ss[i], cfg))
obj <- merge_objects(objs, ss)
obj <- qc_cells(obj, cfg, o$outdir)
obj <- assign_guides(obj, cfg, o$outdir)
obj <- normalize_reduce(obj, cfg, o$outdir)
saveRDS(obj, file.path(o$outdir, "perturbseq_seurat.rds"))
run_dge(obj, cfg, o$outdir)
run_pathways(cfg, o$outdir)
write_report(obj, cfg, o$outdir)
