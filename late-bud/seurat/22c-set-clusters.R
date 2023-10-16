library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

ident.clust <- c("W2.3", "W2.2", "W2.1")
optim.clust <- list("W2.3" = c(5, 4, 3, 2),
                    "W2.2" = c(2),
                    "w2.1" = c(4, 3, 2))

for (i in 1:length(ident.clust)) {
    
    for (t in 1:length(optim.clust[[i]])) {
    
    setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-3")
    
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

path.direct <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-3/", name.direct)
setwd(path.direct)

#filename.robj <- paste0("flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
#save(flox.subclust, file = filename.robj)

flox.subclust.active.ident <- flox.subclust$wildtype.subclust

filename.robj <- paste0("cluster-per-cell-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
save(flox.subclust.active.ident, file = filename.robj)

clust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust$embryo.id, flox.subclust$wildtype.subclust)
colnames(clust.umap) <- c("umap.1", "umap.2", "embryo.id", "wildtype.subclust")

filename.robj <- paste0("umap-by-cluster-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
save(clust.umap, file = filename.robj)

umap.label <- clust.umap %>% group_by(wildtype.subclust) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.subclust), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = wildtype.subclust), size = 3.5) +
labs(title = paste0(ident.clust[i], ", Nipbl Flox/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-cluster-flox-", n, "-clust-", optim.clust[[i]][t], ".png")
ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

name.assay <- paste0("SCT_", n)

marker.genes <- FindAllMarkers(flox.subclust, assay = name.assay, slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

filename.robj <- paste0("marker-genes-flox-", n, "-clust-", optim.clust[[i]][t], ".Robj")
save(marker.genes, file = filename.robj)

filename.csv <- paste0("marker-genes-flox-", n, "-clust-", optim.clust[[i]][t], ".csv")
write.csv(x = marker.genes, file = filename.csv)

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
