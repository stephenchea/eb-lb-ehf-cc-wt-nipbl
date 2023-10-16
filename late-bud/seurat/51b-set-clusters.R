library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

ident.clust <- c("M1", "M2", "M4")
optim.clust <- list("M1" = c(3, 2),
                    "M2" = c(3),
                    "M4" = c(3, 2))

for (i in 1:length(ident.clust)) {
    
    for (t in 1:length(optim.clust[[i]])) {
    
    setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6")
    
filename.robj <- paste0("fin-", ident.clust[i], "-cluster-iters/cluster-per-cell-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
load(filename.robj)

filename.robj <- paste0("fin-", ident.clust[i], "-neighbors.Robj")
load(filename.robj)

renamed.ident <- paste0(ident.clust[i], ".", fin.clust.active.ident)
renamed.ident <- as.character(renamed.ident)
renamed.ident <- factor(renamed.ident)

fin.clust@active.ident <- renamed.ident
Idents(fin.clust) <- renamed.ident
fin.clust$mutant.clust <- fin.clust@active.ident

name.direct <- paste0("fin-", ident.clust[i], "-clust-", optim.clust[[i]][t])
dir.create(name.direct)

path.direct <- paste0("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6/", name.direct)
setwd(path.direct)

fin.clust.active.ident <- fin.clust$mutant.clust

filename.robj <- paste0("cluster-per-cell-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
save(fin.clust.active.ident, file = filename.robj)

filename.robj <- paste0("marker-genes-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
load(file = filename.robj)

top.markers <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/optim.clust[[i]][t])

name.assay <- paste0("SCT_", ident.clust[i])

heatmap <- DoHeatmap(fin.clust, assay = name.assay, slot = "scale.data", features = as.vector(top.markers$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = paste0(ident.clust[i], ", Nipbl FIN/+"), fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("heatmap-marker-genes-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

    }
}
