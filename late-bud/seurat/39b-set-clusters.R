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

load("marker-genes-fin-integr-clust-6.Robj")

marker.genes <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/6)

heatmap <- DoHeatmap(fin.integr, assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(marker.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = "Whole, Nipbl FIN/+", fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "heatmap-marker-genes-fin-integr-clust-6.png", device = "png", units = "in", width = 6.5, height = 9)
