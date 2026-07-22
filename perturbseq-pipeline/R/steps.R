# Author: Mohamed Kassam
# Date: 2026-07-19

import_sample <- function(row, cfg) {
  x <- Read10X_h5(row$h5_path)
  gex_name <- cfg$assays$gene_expression %||% "Gene Expression"
  gex <- if (is.list(x)) x[[gex_name]] else x
  if (is.null(gex)) stop(sprintf("Gene-expression assay '%s' was not found in %s", gex_name, row$h5_path))

  s <- CreateSeuratObject(
    gex,
    project = row$sample_id,
    min.cells = cfg$qc$min_cells_per_gene %||% 3
  )

  # Copy every samplesheet column into cell-level metadata, except the input path.
  metadata_columns <- setdiff(names(row), "h5_path")
  for (column in metadata_columns) {
    s[[column]] <- as.character(row[[column]])
  }

  guide_name <- cfg$assays$guide_capture %||% "CRISPR Guide Capture"
  if (is.list(x) && guide_name %in% names(x)) {
    s[["CRISPR"]] <- CreateAssayObject(x[[guide_name]])
  }
  s
}

merge_objects <- function(objs, ss) {
  if (length(objs) == 1) objs[[1]] else merge(objs[[1]], y = objs[-1], add.cell.ids = ss$sample_id)
}

qc_cells <- function(obj, cfg, out) {
  mito_pattern <- if (identical(tolower(cfg$species %||% "human"), "mouse")) "^mt-" else "^MT-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)
  obj <- CellCycleScoring(
    obj,
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )
  q <- cfg$qc
  obj <- subset(
    obj,
    subset = nFeature_RNA >= q$min_features &
      nFeature_RNA <= q$max_features &
      percent.mt <= q$max_percent_mito
  )
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(file.path(out, "qc_violin.png"), p, width = 12, height = 4)
  obj
}

assign_guides <- function(obj, cfg, out) {
  no_guide <- cfg$guides$no_guide_label %||% "NoGuide"
  if (!"CRISPR" %in% Assays(obj)) {
    obj$guide1 <- no_guide
    obj$guide2 <- NA_character_
    obj$perturbation <- no_guide
    return(obj)
  }

  m <- GetAssayData(obj, assay = "CRISPR", layer = "counts")
  if (nrow(m) == 0) stop("The CRISPR assay contains no guide features")

  top_indices <- apply(m, 2, function(v) order(v, decreasing = TRUE)[seq_len(min(2, length(v)))])
  if (is.null(dim(top_indices))) top_indices <- matrix(top_indices, nrow = 1)

  first_idx <- top_indices[1, ]
  first_umi <- m[cbind(first_idx, seq_len(ncol(m)))]
  first_guide <- rownames(m)[first_idx]

  if (nrow(top_indices) >= 2) {
    second_idx <- top_indices[2, ]
    second_umi <- m[cbind(second_idx, seq_len(ncol(m)))]
    second_guide <- rownames(m)[second_idx]
  } else {
    second_umi <- rep(0, ncol(m))
    second_guide <- rep(NA_character_, ncol(m))
  }

  pass <- first_umi >= cfg$guides$min_umi &
    (second_umi == 0 | first_umi / pmax(second_umi, 1) >= cfg$guides$min_ratio_to_second)

  obj$guide1 <- ifelse(pass, first_guide, no_guide)
  obj$guide2 <- ifelse(pass & second_umi >= cfg$guides$min_umi, second_guide, NA_character_)
  target <- function(x) sub("([_-]g?[0-9]+)$", "", x)
  obj$perturbation <- ifelse(pass, target(first_guide), no_guide)
  fwrite(as.data.frame(obj[[]]), file.path(out, "cell_metadata.tsv"), sep = "\t")
  obj
}

normalize_reduce <- function(obj, cfg, out) {
  regress <- intersect(unlist(cfg$normalization$regress), colnames(obj[[]]))
  obj <- SCTransform(
    obj,
    vst.flavor = "v2",
    vars.to.regress = regress,
    variable.features.n = cfg$normalization$variable_features,
    verbose = FALSE
  )
  dims <- seq_len(cfg$reduction$npcs)
  obj <- RunPCA(obj, npcs = cfg$reduction$npcs, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cfg$reduction$resolution, verbose = FALSE)

  configured_groups <- unlist(cfg$visualization$group_by %||% c("treatment", "perturbation"))
  group_by <- intersect(configured_groups, colnames(obj[[]]))
  if (length(group_by) > 0) {
    p <- DimPlot(obj, reduction = "umap", group.by = group_by)
    ggsave(file.path(out, "umap_conditions.png"), p, width = 12, height = 5)
  }
  obj
}

run_dge <- function(obj, cfg, out) {
  d <- safe_dir(file.path(out, "dge"))
  identity_column <- cfg$dge$identity_column %||% "perturbation"
  if (!identity_column %in% colnames(obj[[]])) {
    stop(sprintf("Configured DGE identity column '%s' is missing from cell metadata", identity_column))
  }
  Idents(obj) <- identity_column
  levels_present <- levels(Idents(obj))
  control <- cfg$dge$reference_level %||% cfg$references$perturbation %||% cfg$guides$no_guide_label
  if (is.null(control) || !control %in% levels_present) {
    stop(sprintf("Configured DGE reference level '%s' is not present. Available levels: %s", control, paste(levels_present, collapse = ", ")))
  }

  for (level in setdiff(levels_present, control)) {
    if (sum(Idents(obj) == level) < (cfg$dge$min_cells_per_group %||% 10)) next
    z <- FindMarkers(
      obj,
      ident.1 = level,
      ident.2 = control,
      assay = "SCT",
      test.use = cfg$dge$method,
      min.pct = cfg$dge$min_pct,
      logfc.threshold = cfg$dge$logfc_threshold
    )
    z$gene <- rownames(z)
    z$contrast <- paste0(level, "_vs_", control)
    fwrite(z, file.path(d, paste0(make.names(z$contrast[1]), ".tsv")), sep = "\t")
  }
}

run_pathways <- function(cfg, out) {
  if (!isTRUE(cfg$pathway$enabled)) return(invisible(NULL))
  suppressPackageStartupMessages({
    library(fgsea)
    library(msigdbr)
  })

  d <- file.path(out, "dge")
  fs <- list.files(d, "\\.tsv$", full.names = TRUE)
  if (!length(fs)) return(invisible(NULL))
  pathway_dir <- safe_dir(file.path(out, "pathway"))

  for (f in fs) {
    z <- fread(f)
    lfc_candidates <- intersect(c("avg_log2FC", "logFC"), names(z))
    if (!length(lfc_candidates) || !"gene" %in% names(z)) next
    ranks <- setNames(z[[lfc_candidates[1]]], z$gene)
    ranks <- sort(ranks[is.finite(ranks) & !is.na(names(ranks))], decreasing = TRUE)
    if (!length(ranks)) next

    species_name <- if (identical(tolower(cfg$species %||% "human"), "mouse")) "Mus musculus" else "Homo sapiens"
    collections <- unlist(cfg$pathway$collections %||% "H")
    for (collection in collections) {
      ms <- msigdbr(species = species_name, category = collection)
      pathways <- split(ms$gene_symbol, ms$gs_name)
      fg <- fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = cfg$pathway$min_size,
        maxSize = cfg$pathway$max_size
      )
      output_name <- paste0(tools::file_path_sans_ext(basename(f)), "_", make.names(collection), ".tsv")
      fwrite(fg, file.path(pathway_dir, output_name), sep = "\t")
    }
  }
}

write_report <- function(obj, cfg, out) {
  image_html <- if (file.exists(file.path(out, "umap_conditions.png"))) '<img src="umap_conditions.png" width="1000">' else ""
  html <- paste0(
    "<html><body><h1>Perturb-seq Analysis Report</h1>",
    "<p>Cells retained: ", ncol(obj), "</p>",
    "<p>Genes: ", nrow(obj), "</p>",
    "<p>Generated: ", Sys.time(), "</p>",
    image_html,
    "</body></html>"
  )
  writeLines(html, file.path(out, "report.html"))
}
