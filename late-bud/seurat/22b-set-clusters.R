library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

ident.clust <- c("W1.1")
optim.clust <- list("W1.1" = c(3))

for (i in 1:length(ident.clust)) {
    
    for (t in 1:length(optim.clust[[i]])) {
    
    setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W1-clust-2")
    
    n <- str_replace_all(ident.clust[i], "[.]", "_")
    
filename.robj <- paste0("flox-", n, "-cluster-iters/cluster-per-cell-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
load(filename.robj)

filename.robj <- paste0("flox-", n, "-neighbors.Robj")
load(filename.robj)

renamed.ident <- paste0(ident.clust[i], ".", flox.subclust.active.ident)
renamed.ident <- as.character(renamed.ident)
renamed.ident <- factor(renamed.ident)

flox.subclust@active.ident <- renamed.ident
Idents(flox.subclust) <- renamed.ident
flox.subclust$wildtype.subclust <- flox.subclust@active.ident

name.direct <- paste0("flox-", n, "-clust-", optim.clust[[i]][t])
dir.create(name.direct)

path.direct <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W1-clust-2/", name.direct)
setwd(path.direct)

flox.subclust.active.ident <- flox.subclust$wildtype.subclust

filename.robj <- paste0("cluster-per-cell-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
save(flox.subclust.active.ident, file = filename.robj)

name.assay <- paste0("SCT_", n)

filename.robj <- paste0("marker-genes-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
load(file = filename.robj)

top.markers <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/optim.clust[[i]][t])

heatmap <- DoHeatmap(flox.subclust, assay = name.assay, slot = "scale.data", features = as.vector(top.markers$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = paste0(ident.clust[i], ", Nipbl Flox/+"), fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("heatmap-marker-genes-flox-", n, "-clust-", optim.clust[[i]][t], ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

    }
}


ident.clust <- c("W2.1", "W2.2")
optim.clust <- list("W2.1" = c(2),
                    "W2.2" = c(2))

for (i in 1:length(ident.clust)) {
    
    for (t in 1:length(optim.clust[[i]])) {
    
    setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-2")
    
    n <- str_replace_all(ident.clust[i], "[.]", "_")
    
filename.robj <- paste0("flox-", n, "-cluster-iters/cluster-per-cell-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
load(filename.robj)

filename.robj <- paste0("flox-", n, "-neighbors.Robj")
load(filename.robj)

renamed.ident <- paste0(ident.clust[i], ".", flox.subclust.active.ident)
renamed.ident <- as.character(renamed.ident)
renamed.ident <- factor(renamed.ident)

flox.subclust@active.ident <- renamed.ident
Idents(flox.subclust) <- renamed.ident
flox.subclust$wildtype.subclust <- flox.subclust@active.ident

name.direct <- paste0("flox-", n, "-clust-", optim.clust[[i]][t])
dir.create(name.direct)

path.direct <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-2/", name.direct)
setwd(path.direct)

flox.subclust.active.ident <- flox.subclust$wildtype.subclust

filename.robj <- paste0("cluster-per-cell-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
save(flox.subclust.active.ident, file = filename.robj)

name.assay <- paste0("SCT_", n)

filename.robj <- paste0("marker-genes-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
load(file = filename.robj)

top.markers <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/optim.clust[[i]][t])

heatmap <- DoHeatmap(flox.subclust, assay = name.assay, slot = "scale.data", features = as.vector(top.markers$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = paste0(ident.clust[i], ", Nipbl Flox/+"), fill = "Norm\nUMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("heatmap-marker-genes-flox-", n, "-clust-", optim.clust[[i]][t], ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

    }
}
