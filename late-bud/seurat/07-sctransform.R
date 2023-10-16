print("load libraries")
library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat")

print("load seurat object")
load("all-aggr-seurat.Robj")

all.aggr$experiment <- factor(all.aggr$experiment, levels = c("20171003", "20180515", "20181220"))
all.aggr$embryo.id <- factor(all.aggr$embryo.id, levels = c("AVJ4-56", "AYQ8-11", "AYR4-57", "BBH6-17", "AXL7-21", "ATU6-45", "ATU6-69", "AYR4-22", "AZK2-45", "BBH6-59", "AXL7-69"))
all.aggr$genotype <- factor(all.aggr$genotype, levels = c("Flox", "FIN"))

print("load outlier cells")
load("all-out-all-aggr.Robj")

print("add outlier cells into seurat object")
all.aggr$all.out <- all.out

print("set identities to outlier")
Idents(all.aggr) <- all.aggr$all.out

print("subset high mito cells")
all.aggr <- subset(all.aggr, idents = "FALSE")
dim(all.aggr)

print("apply sctransform")
all.aggr <- SCTransform(object = all.aggr, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 1, n_cells = round(dim(all.aggr)[2]*.5), n_genes = round(dim(all.aggr)[1]*.1), new.assay.name = "SCT_ALL", verbose = TRUE)

save(object = all.aggr, file = "all-aggr-sctrans.Robj")
