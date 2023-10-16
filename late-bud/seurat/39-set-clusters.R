library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat")
load("fin-integr-neighbors.Robj")

########################################
#set number of clusters to 2

load("fin-integr-cluster-iters/cluster-per-cell-fin-integr-clust-6.Robj")

fin.integr@active.ident <- fin.integr.active.ident
fin.integr <- RenameIdents(object = fin.integr, "1" = "M1", "2" = "M2", "3" = "M3", "4" = "M4", "5" = "M5", "6" = "M6")
fin.integr@active.ident <- factor(fin.integr@active.ident, levels = c("M1", "M2", "M3", "M4", "M5", "M6"))
fin.integr$mutant.whole <- fin.integr@active.ident

dir.create("fin-integr-clust-6")
setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6")

#save(fin.integr, file = "fin-integr-clust-6.Robj")

fin.integr.active.ident <- fin.integr$mutant.whole
save(fin.integr.active.ident, file = "cluster-per-cell-fin-integr-clust-6.Robj")

whole.umap <- data.frame(fin.integr@reductions$umap@cell.embeddings, fin.integr$embryo.id, fin.integr$mutant.whole)
colnames(whole.umap) <- c("umap.1", "umap.2", "embryo.id", "mutant.whole")

save(whole.umap, file = "umap-by-cluster-fin-integr-clust-6.Robj")

umap.label <- whole.umap %>% group_by(mutant.whole) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = mutant.whole), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = mutant.whole), size = 3.5) +
labs(title = "Whole, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "umap-by-cluster-fin-integr-clust-6.png", device = "png", width = 6.5, height = 4.5, units = "in")

marker.genes <- FindAllMarkers(fin.integr, assay = "SCT_GENOTYPE", slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

save(marker.genes, file = "marker-genes-fin-integr-clust-6.Robj")

write.csv(x = marker.genes, file = "marker-genes-fin-integr-clust-6.csv")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice(1:10)

heatmap <- DoHeatmap(fin.integr, assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "Whole, Nipbl FIN/+", x = "Cluster", y = "Top Genes by AUC", fill = "Centered\nNorm UMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-fin-integr-clust-6.png", device = "png", units = "in", width = 6.5, height = 9)
