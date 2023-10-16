library(Seurat)
library(tidyverse)
library(RColorBrewer)

setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18")

load("fin-integr-projected-flox-integr-clust-5-18.Robj")

fin.integr.active.ident <- fin.integr$mutant.whole
save(fin.integr.active.ident, file = "cluster-per-cell-by-mutant-fin-integr-projected-flox-integr-clust-5-18.Robj")

fin.integr.active.ident <- fin.integr$wildtype.whole
save(fin.integr.active.ident, file = "cluster-per-cell-by-wildtype-fin-integr-projected-flox-integr-clust-5-18.Robj")

whole.umap <- data.frame(fin.integr@reductions$ref.umap@cell.embeddings, fin.integr$embryo.id, fin.integr$mutant.whole, fin.integr$wildtype.whole)
colnames(whole.umap) <- c("umap.1", "umap.2", "embryo.id", "mutant.whole", "wildtype.whole")

save(whole.umap, file = "projected-umap-by-cluster-fin-integr-projected-flox-integr-clust-5-18.Robj")

whole.label <- whole.umap %>% group_by(mutant.whole) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = mutant.whole), size = 0.1) +
geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = mutant.whole), size = 3.5) +
labs(title = "Whole, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "projected-umap-by-mutant-fin-integr-projected-flox-integr-clust-5-18.png", device = "png", width = 6.5, height = 4.5, units = "in")

whole.label <- whole.umap %>% group_by(wildtype.whole) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.whole), size = 0.1) +
geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = wildtype.whole), size = 3.5) +
labs(title = "Whole, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "projected-umap-by-wildtype-fin-integr-projected-flox-integr-clust-5-18.png", device = "png", width = 6.5, height = 4.5, units = "in")

Idents(fin.integr) <- fin.integr$wildtype.whole

fin.integr <- PrepSCTFindMarkers(fin.integr, assay = "SCT_GENOTYPE", verbose = TRUE)

marker.genes <- FindAllMarkers(fin.integr, assay = "SCT_GENOTYPE", slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

save(marker.genes, file = "marker-genes-fin-integr-projected-flox-integr-clust-5-18.Robj")

write.csv(x = marker.genes, file = "marker-genes-fin-integr-projected-flox-integr-clust-5-18.csv")

top.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/15)

heatmap <- DoHeatmap(fin.integr, assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = FALSE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "Whole, Nipbl FIN/+", fill = "Norm\nUMI", color = "Cluster") +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-fin-integr-projected-flox-integr-clust-5-18.png", device = "png", units = "in", width = 6.5, height = 9)
