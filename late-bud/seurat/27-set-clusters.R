library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-2/flox-W2_1-clust-2")

load("flox-W2_1-clust-2.Robj")
load("umap-by-cluster-flox-W2_1-clust-2-5.Robj")

flox.subclust$wildtype.subclust <- subclust.umap$wildtype.subclust
Idents(flox.subclust) <- subclust.umap$wildtype.subclust

flox.subclust.active.ident <- flox.subclust$wildtype.subclust
save(flox.subclust.active.ident, file = "cluster-per-cell-flox-W2_1-clust-2-5.Robj")

subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust$embryo.id, flox.subclust$wildtype.subclust)
colnames(subclust.umap) <- c("umap.1", "umap.2", "embryo.id", "wildtype.subclust")

save(subclust.umap, file = "umap-by-cluster-flox-W2_1-clust-2-5.Robj")

umap.label <- subclust.umap %>% group_by(wildtype.subclust) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.subclust), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = wildtype.subclust), size = 3.5) +
labs(title = "W2_1, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "umap-by-cluster-flox-W2_1-clust-2-5.png", device = "png", width = 6.5, height = 4.5, units = "in")

marker.genes <- FindAllMarkers(flox.subclust, assay = "SCT_W2_1", slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

save(marker.genes, file = "marker-genes-flox-W2_1-clust-2-5.Robj")

write.csv(x = marker.genes, file = "marker-genes-flox-W2_1-clust-2-5.csv")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/5)

heatmap <- DoHeatmap(flox.subclust, assay = "SCT_W2_1", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "W2.1, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-flox-W2_1-clust-2-5.png", device = "png", units = "in", width = 6.5, height = 9)
