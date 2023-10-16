library(Seurat)
library(tidyverse)
library(RColorBrewer)

load("/dfs6/pub/schea2/20210505-seurat/flox-integr-neighbors.Robj")

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5")
load("umap-by-germ-layer-flox-integr-clust-5-18.Robj")

whole.umap <- whole.umap
flox.integr$germ.whole <- whole.umap$germ.whole
Idents(flox.integr) <- flox.integr$germ.whole

flox.integr <- PrepSCTFindMarkers(flox.integr, assay = "SCT_GENOTYPE", verbose = TRUE)

diff.genes <- FindAllMarkers(flox.integr, assay = "SCT_GENOTYPE", slot = "data", only.pos = TRUE, test.use = "wilcox", min.pct = 0, min.diff.pct = -Inf, logfc.threshold = 0, return.thresh = 0.5, max.cells.per.ident = dim(flox.integr)[2]/length(levels(flox.integr$germ.whole)), verbose = TRUE)

save(diff.genes, file = "diff-genes-by-germ-layer-flox-integr-clust-5-18.Robj")

write.csv(x = diff.genes, file = "diff-genes-by-germ-layer-flox-integr-clust-5-18.csv")

top.genes <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(flox.integr$germ.whole)))

heatmap <- DoHeatmap(flox.integr, assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "Whole, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-diff-genes-by-germ-layer-flox-integr-clust-5-18.png", device = "png", units = "in", width = 6.5, height = 9)
