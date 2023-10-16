print("loading libraries")
library(Seurat)
library(tidyverse)

print("setting working directory")
setwd("/dfs6/pub/schea2/20210505-seurat")

print("load seurat object")
load("all-aggr-sctrans.Robj")

print("perform principal component analysis")
all.aggr <- RunPCA(object = all.aggr, assay = "SCT_ALL", features = all.aggr@assays$SCT_ALL@var.features, npcs = 50, verbose = TRUE)

print("save seurat object")
save(object = all.aggr, file = "all-aggr-pca.Robj")

print("save standard deviation of principal components")
stdev.pca <- all.aggr@reductions$pca@stdev
save(object = stdev.pca, file = "stdev-pc-all-aggr.Robj")
