library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6/fin-M2-clust-3")

load("fin-M2-clust-3.Robj")
load("umap-by-cluster-fin-M2-clust-3-8.Robj")

fin.clust$mutant.clust <- clust.umap$mutant.clust
Idents(fin.clust) <- clust.umap$mutant.clust

fin.clust.active.ident <- fin.clust$mutant.clust
save(fin.clust.active.ident, file = "cluster-per-cell-fin-M2-clust-3-8.Robj")

clust.umap <- data.frame(fin.clust@reductions$umap@cell.embeddings, fin.clust$embryo.id, fin.clust$mutant.clust)
colnames(clust.umap) <- c("umap.1", "umap.2", "embryo.id", "mutant.clust")

save(clust.umap, file = "umap-by-cluster-fin-M2-clust-3-8.Robj")

umap.label <- clust.umap %>% group_by(mutant.clust) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = mutant.clust), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = mutant.clust), size = 3.5) +
labs(title = "M2, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "umap-by-cluster-fin-M2-clust-3-8.png", device = "png", width = 6.5, height = 4.5, units = "in")

marker.genes <- FindAllMarkers(fin.clust, assay = "SCT_M2", slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

save(marker.genes, file = "marker-genes-fin-M2-clust-3-8.Robj")

write.csv(x = marker.genes, file = "marker-genes-fin-M2-clust-3-8.csv")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/8)

heatmap <- DoHeatmap(fin.clust, assay = "SCT_M2", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "M2, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-fin-M2-clust-3-8.png", device = "png", units = "in", width = 6.5, height = 9)
