library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)

setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

load("/share/crsp/lab/alcalof/schea2/20220414-seurat/all-aggr-neighbors.Robj")
Idents(all.aggr) <- all.aggr$genotype

#subset flox cells from all cells
flox.aggr <- subset(all.aggr, idents = "Nipbl Flox/+")

remove(all.aggr)

flox.aggr <- SCTransform(object = flox.aggr, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 1, n_cells = round(mean(table(flox.aggr$label.embryo))), n_genes = round(mean(flox.aggr$nFeature_RNA)), new.assay.name = "SCT_GENOTYPE", verbose = TRUE)

Idents(flox.aggr) <- flox.aggr$label.embryo

#split flox cells according to sample id
embryo.flox <- SplitObject(flox.aggr, split.by = "label.embryo")

remove(flox.aggr)

for (i in 1:length(embryo.flox)) {
    
    embryo.flox[[i]] <- SCTransform(embryo.flox[[i]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = TRUE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.flox[[i]]$nFeature_RNA)), verbose = TRUE)
    
    DefaultAssay(embryo.flox[[i]]) <- "SCT_EMBRYO"
    
    }

save(embryo.flox, file = "flox-aggr-sctrans.Robj")

print("grabbing variable genes from each embryo")
var.genes <- list(embryo.flox[[1]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[2]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[3]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[4]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[5]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[6]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[7]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[8]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[9]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[10]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[11]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[12]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[13]]@assays$SCT_EMBRYO@var.features)

names(var.genes) <- names(embryo.flox)

save(var.genes, file = "var-genes-flox-aggr-sctrans.Robj")
