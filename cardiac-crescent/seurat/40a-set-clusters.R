library(Seurat)
library(tidyverse)

setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-neighbors.Robj")

res <- as.character(c("0.0259", "0.0151", "0.0092", "0.0065"))

for (i in 1:length(res)) {

filename.robj <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-cluster-iters/cluster-per-cell-fin-integr-res-", res[i], ".Robj")

load(filename.robj)

fin.integr@active.ident <- fin.integr.active.ident
fin.integr$cluster.whole <- fin.integr@active.ident

path.direct <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-", res[i])
dir.create(path.direct)
setwd(path.direct)

cluster.whole <- fin.integr$cluster.whole
filename.robj <- paste0("cluster-per-cell-fin-integr-res-", res[i], ".Robj")
save(cluster.whole, file = filename.robj)

umap.whole <- data.frame(fin.integr@reductions$umap@cell.embeddings, fin.integr$stage, fin.integr$label.embryo, fin.integr$cluster.whole)
colnames(umap.whole) <- c("umap.1", "umap.2", "stage", "embryo", "cluster.whole")

filename.robj <- paste0("umap-by-cluster-fin-integr-res-", res[i], ".Robj")
save(umap.whole, file = filename.robj)

#label.whole <- umap.whole %>% group_by(cluster.whole) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

#cluster.umap <- ggplot() +
#geom_point(data = sample(umap.whole), mapping = aes(x = umap.1, y = umap.2, color = cluster.whole), shape = 20, stroke = 0, size = 0.5) +
#geom_text(data = label.whole, mapping = aes(x = umap.1, y = umap.2, label = cluster.whole), size = 2.3) +
#labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
#guides(color = guide_legend(override.aes = list(size = 2.3))) +
#theme_classic(base_size = 7) +
#theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

#filename.png <- paste0("umap-by-cluster-fin-integr-res-", res[i], ".png")
#ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

diff.genes <- FindAllMarkers(fin.integr, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.integr)[2]/length(levels(fin.integr$cluster.whole)), return.thresh = 0.05)

filename.robj <- paste0("diff-genes-fin-integr-res-", res[i], ".Robj")
save(diff.genes, file = filename.robj)

filename.csv <- paste0("diff-genes-fin-integr-res-", res[i], ".csv")
write.csv(x = diff.genes, file = filename.csv)

#top.markers <- diff.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 177/length(levels(fin.integr$cluster.whole)))

#cell.names <- rownames(fin.integr@meta.data)

#heatmap <- DoHeatmap(fin.integr, cells = sample(cell.names, round(length(cell.names)*0.5), replace = FALSE), assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(top.markers$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
#scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
#scale_x_discrete(position = "top") +
#labs(title = "Nipbl Flox/+", x = "Cluster", y = "Top Genes", fill = "Standardized\nLog1p UMI") +
#guides(color = FALSE) +
#theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

#filename.png <- paste0("heatmap-diff-genes-fin-integr-res-", res[i], ".png")
#ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 9, height = 23.5)

setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

}
