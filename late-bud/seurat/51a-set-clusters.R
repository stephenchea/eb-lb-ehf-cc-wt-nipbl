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

#filename.robj <- paste0("fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
#save(fin.clust, file = filename.robj)

fin.clust.active.ident <- fin.clust$mutant.clust

filename.robj <- paste0("cluster-per-cell-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
save(fin.clust.active.ident, file = filename.robj)

clust.umap <- data.frame(fin.clust@reductions$umap@cell.embeddings, fin.clust$embryo.id, fin.clust$mutant.clust)
colnames(clust.umap) <- c("umap.1", "umap.2", "embryo.id", "mutant.clust")

filename.robj <- paste0("umap-by-cluster-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
save(clust.umap, file = filename.robj)

filename.csv <- paste0("umap-by-cluster-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".csv")
write.csv(x = clust.umap, file = filename.csv)

umap.label <- clust.umap %>% group_by(mutant.clust) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = mutant.clust), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = mutant.clust), size = 3.5) +
labs(title = paste0(ident.clust[i], ", Nipbl FIN/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-cluster-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".png")
ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

name.assay <- paste0("SCT_", ident.clust[i])

marker.genes <- FindAllMarkers(fin.clust, assay = name.assay, slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

filename.robj <- paste0("marker-genes-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".Robj")
save(marker.genes, file = filename.robj)

filename.csv <- paste0("marker-genes-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".csv")
write.csv(x = marker.genes, file = filename.csv)

top.markers <- marker.genes %>% group_by(cluster) %>% slice(1:10)

heatmap <- DoHeatmap(fin.clust, assay = name.assay, slot = "scale.data", features = as.vector(top.markers$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = paste0(ident.clust[i], ", Nipbl FIN/+"), x = "Cluster", y = "Top Genes by AUC", fill = "Centered\nNorm UMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("heatmap-marker-genes-fin-", ident.clust[i], "-clust-", optim.clust[[i]][t], ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

    }
}
