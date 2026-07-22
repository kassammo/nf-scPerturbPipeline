# Author: Mohamed Kassam
# Date: 2026-07-19
suppressPackageStartupMessages({library(Seurat);library(data.table);library(ggplot2)})
`%||%` <- function(x,y) if(is.null(x)) y else x
safe_dir <- function(x){dir.create(x,recursive=TRUE,showWarnings=FALSE);x}
write_tsv <- function(x,path) fwrite(as.data.frame(x),path,sep='\t')
