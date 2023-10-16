library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat")
load("flox-integr-neighbors.Robj")

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5")
load("cluster-per-cell-flox-integr-clust-5.Robj")

flox.integr$wildtype.whole <- flox.integr.active.ident
Idents(flox.integr) <- flox.integr.active.ident

load("marker-genes-flox-integr-clust-5.Robj")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/length(levels(flox.integr$wildtype.whole)))

heatmap <- DoHeatmap(flox.integr, assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "Whole, Nipbl Flox/+", fill = "Norm\nUMI") +
guides(color = "none") +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-flox-integr-clust-5.png", device = "png", units = "in", width = 6.5, height = 9)
