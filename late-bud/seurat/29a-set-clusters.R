library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W1-clust-2")

load("flox-W1-clust-2.Robj")
load("umap-by-cluster-flox-W1-clust-2-4.Robj")

flox.clust$wildtype.clust <- clust.umap$wildtype.clust
Idents(flox.clust) <- clust.umap$wildtype.clust

flox.clust.active.ident <- flox.clust$wildtype.clust
save(flox.clust.active.ident, file = "cluster-per-cell-flox-W1-clust-2-4.Robj")

clust.umap <- data.frame(flox.clust@reductions$umap@cell.embeddings, flox.clust$embryo.id, flox.clust$wildtype.clust)
colnames(clust.umap) <- c("umap.1", "umap.2", "embryo.id", "wildtype.clust")

save(clust.umap, file = "umap-by-cluster-flox-W1-clust-2-4.Robj")

umap.label <- clust.umap %>% group_by(wildtype.clust) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.clust), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = wildtype.clust), size = 3.5) +
labs(title = "W1, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "umap-by-cluster-flox-W1-clust-2-4.png", device = "png", width = 6.5, height = 4.5, units = "in")

marker.genes <- FindAllMarkers(flox.clust, assay = "SCT_GENOTYPE", slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

save(marker.genes, file = "marker-genes-flox-W1-clust-2-4.Robj")

write.csv(x = marker.genes, file = "marker-genes-flox-W1-clust-2-4.csv")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/4)

heatmap <- DoHeatmap(flox.clust, assay = "SCT_W1", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "W1, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-flox-W1-clust-2-4.png", device = "png", units = "in", width = 6.5, height = 9)

########################################
setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-2")

load("flox-W2-clust-2.Robj")
load("umap-by-cluster-flox-W2-clust-2-6.Robj")

flox.clust$wildtype.clust <- umap.clust$wildtype.clust
Idents(flox.clust) <- umap.clust$wildtype.clust

flox.clust.active.ident <- flox.clust$wildtype.clust
save(flox.clust.active.ident, file = "cluster-per-cell-flox-W2-clust-2-6.Robj")

clust.umap <- data.frame(flox.clust@reductions$umap@cell.embeddings, flox.clust$embryo.id, flox.clust$wildtype.clust)
colnames(clust.umap) <- c("umap.1", "umap.2", "embryo.id", "wildtype.clust")

save(clust.umap, file = "umap-by-cluster-flox-W2-clust-2-6.Robj")

umap.label <- clust.umap %>% group_by(wildtype.clust) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.clust), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = wildtype.clust), size = 3.5) +
labs(title = "W2, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "umap-by-cluster-flox-W2-clust-2-6.png", device = "png", width = 6.5, height = 4.5, units = "in")

marker.genes <- FindAllMarkers(flox.clust, assay = "SCT_W2", slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

save(marker.genes, file = "marker-genes-flox-W2-clust-2-6.Robj")

write.csv(x = marker.genes, file = "marker-genes-flox-W2-clust-2-6.csv")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/6)

heatmap <- DoHeatmap(flox.clust, assay = "SCT_W2", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "W2, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-flox-W2-clust-2-6.png", device = "png", units = "in", width = 6.5, height = 9)
