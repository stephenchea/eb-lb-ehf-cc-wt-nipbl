print("load libraries")
library(Seurat)
library(pracma)

setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

load("all-aggr-seurat.Robj")
load("all-out-all-aggr.Robj")

print("add outlier cells into seurat object")
all.aggr$all.out <- all.out

print("set identities to outlier")
Idents(all.aggr) <- all.aggr$all.out

print("subset high quality cells")
all.aggr <- subset(all.aggr, idents = "FALSE")
dim(all.aggr)

print("apply sctransform")
all.aggr <- SCTransform(object = all.aggr, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 1, n_cells = mean(table(all.aggr$label.embryo)), n_genes = mean(all.aggr$nFeature_RNA), new.assay.name = "SCT_ALL", verbose = TRUE)

all.aggr <- RunPCA(object = all.aggr, assay = "SCT_ALL", features = all.aggr@assays$SCT_ALL@var.features, npcs = 50, verbose = TRUE)

stdev.pca <- all.aggr@reductions$pca@stdev
dens.stdev <- density(stdev.pca)
peak.stdev <- findpeaks(x = dens.stdev$y*-1, npeaks = 1)[1, 2]
elbow.stdev <- dens.stdev$x[peak.stdev]
pc.number <- sum(stdev.pca > elbow.stdev)
print(pc.number)

save(pc.number, file = "sig-pc-all-aggr.Robj")

all.aggr <- RunUMAP(all.aggr, reduction = "pca", dims = 1:pc.number)

whole.umap <- data.frame(all.aggr@reductions$umap@cell.embeddings, all.aggr$stage, all.aggr$label.embryo, all.aggr$genotype)
colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "embryo", "genotype")
save(whole.umap, file = "umap-by-embryo-all-aggr.Robj")

all.aggr <- FindNeighbors(all.aggr, dims = 1:pc.number)

save(all.aggr, file = "all-aggr-neighbors.Robj")
