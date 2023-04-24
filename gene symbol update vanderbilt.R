
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(AnnotationHub)
library(ensembldb)
library(cowplot)
library(ggplot2)
library(scales)
library(Seurat.utils)

counts <- readMM("matrix.mtx")

# Read in `genes.tsv`
genes <- read_tsv("genes.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cells <- read_tsv("GSE135893_barcodes.tsv.gz", col_names = FALSE)
cell_ids=cells$X1
counts <- as(counts, "dgCMatrix")
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

meta=read.csv("metadata.csv", header=T, row.names = 1)

vand=CreateSeuratObject(counts = counts)
vand=UpdateGenesSeurat(vand)

vand@meta.data=meta
DotPlot(vand, features="RACK1", group.by="orig.ident")
DotPlot(vand, features="GNB2L1")
saveRDS(vand, file="vand.rds")
