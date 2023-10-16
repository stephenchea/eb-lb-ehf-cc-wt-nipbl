print("loading libraries")
library(Seurat)
library(tidyverse)

print("loading neighbors")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/all-aggr-neighbors.Robj")

print("creating cluster iters folder")
dir.create("all-aggr-cluster-iters")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/all-aggr-cluster-iters")

for (i in seq(from = 0, to = 0.5, by = 0.0001)) {

print("res")
print(i)

print("finding clusters")
all.aggr <- FindClusters(all.aggr, resolution = i, verbose = TRUE)

print("renumbering clusters")
renumbered.cluster.ids <- 1:length(levels(all.aggr))
names(renumbered.cluster.ids) <- levels(all.aggr)
all.aggr <- RenameIdents(all.aggr, renumbered.cluster.ids)
all.aggr.active.ident <- all.aggr@active.ident

print("saving cluster identities")
filename.robj <- paste0("cluster-per-cell-all-aggr-res-", sprintf(fmt = "%.4f", i), ".Robj")
save(all.aggr.active.ident, file = filename.robj)

print("generating umap coordinates")
whole.umap <- data.frame(all.aggr@reductions$umap@cell.embeddings, all.aggr$stage, all.aggr$genotype, all.aggr$label.embryo, all.aggr@active.ident)
colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "genotype", "embryo", "cluster")

print("generating labels for umap")
whole.label <- whole.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

print("generating umap plot")
cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), shape = 20, stroke = 0, size = 0.5) +
geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 2.3) +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

print("saving umap plot")
filename.png <- paste0("umap-by-cluster-all-aggr-res-", sprintf(fmt = "%.4f", i), ".png")
#ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

}

