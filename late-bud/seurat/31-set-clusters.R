library(Seurat)
library(tidyverse)
library(RColorBrewer)

load("/dfs6/pub/schea2/20210505-seurat/flox-integr-neighbors.Robj")

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5")
load("umap-by-cluster-flox-integr-clust-5-18.Robj")

whole.umap <- whole.umap
flox.integr$wildtype.whole <- whole.umap$wildtype.whole
Idents(flox.integr) <- flox.integr$wildtype.whole

flox.integr.active.ident <- flox.integr$wildtype.whole
save(flox.integr.active.ident, file = "cluster-per-cell-flox-integr-clust-5-18.Robj")

whole.umap <- data.frame(flox.integr@reductions$umap@cell.embeddings, flox.integr$embryo.id, flox.integr$wildtype.whole)
colnames(whole.umap) <- c("umap.1", "umap.2", "embryo.id", "wildtype.whole")

save(whole.umap, file = "umap-by-cluster-flox-integr-clust-5-18.Robj")

umap.label <- whole.umap %>% group_by(wildtype.whole) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.whole), size = 0.1) +
geom_text(data = umap.label %>% filter, mapping = aes(x = umap.1, y = umap.2, label = wildtype.whole), size = 3.5) +
labs(title = "Whole, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "umap-by-cluster-flox-integr-clust-5-18.png", device = "png", width = 6.5, height = 4.5, units = "in")

flox.integr <- PrepSCTFindMarkers(flox.integr, assay = "SCT_GENOTYPE", verbose = TRUE)

diff.genes <- FindAllMarkers(flox.integr, assay = "SCT_GENOTYPE", slot = "data", only.pos = TRUE, test.use = "wilcox", min.pct = 0, min.diff.pct = -Inf, logfc.threshold = 0, return.thresh = 0.5, max.cells.per.ident = dim(flox.integr)[2]/length(levels(flox.integr$wildtype.whole)), verbose = TRUE)

save(diff.genes, file = "diff-genes-flox-integr-clust-5-18.Robj")

write.csv(x = diff.genes, file = "diff-genes-flox-integr-clust-5-18.csv")

top.genes <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(flox.integr$wildtype.whole)))

heatmap <- DoHeatmap(flox.integr, assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "Whole, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-flox-integr-clust-5-18.png", device = "png", units = "in", width = 6.5, height = 9)
