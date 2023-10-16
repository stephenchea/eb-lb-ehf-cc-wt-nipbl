print("loading libraries")
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(patchwork)
setwd("/dfs6/pub/schea2/20210505-seurat")

print("perform principal component analysis")
load("all-aggr-pca.Robj")

print("save standard deviation of principal components")
stdev.pca <- all.aggr@reductions$pca@stdev
save(object = stdev.pca, file = "stdev-pc-all-aggr.Robj")

var.exp <- stdev.pca^2/sum(stdev.pca^2)

print("set number of principal components")
d.var.exp <- var.exp - median(var.exp)
pc.number <- sum(d.var.exp > mad(var.exp)*2)
save(pc.number, file = "sig-pc-all-aggr.Robj")

all.aggr <- RunTSNE(all.aggr, dims = 1:pc.number)

tsne.whole <- data.frame(all.aggr@reductions$tsne@cell.embeddings, all.aggr$embryo.id, all.aggr$genotype)
colnames(tsne.whole) <- c("tsne.1", "tsne.2", "embryo.id", "genotype")

tsne.whole$genotype <- factor(tsne.whole$genotype, levels = c("Flox", "FIN"))

save(tsne.whole, file = "tsne-coord-all-aggr-tsne.Robj")

tsne.genotype <- ggplot() +
geom_point(data = sample(tsne.whole), mapping = aes(x = tsne.1, y = tsne.2, color = genotype), size = 0.1) +
scale_color_manual(values = c("royalblue", "indianred"), labels = c("Nipbl Flox/+", "Nipbl FIN/+")) +
labs(title = "Nipbl Flox/+ & FIN/+", x = "t-SNE 1", y = "t-SNE 2", color = "Genotype") +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = tsne.genotype, filename = "tsne-by-genotype-all-aggr.png", device = "png", width = 6.5, height = 4.5)

tsne.flox <- ggplot() +
geom_point(data = sample(tsne.whole %>% filter(genotype == "Flox")), mapping = aes(x = tsne.1, y = tsne.2, color = embryo.id), size = 0.1) +
scale_color_manual(labels = c("EP1", "EP2", "EP3", "EP4", "EP5"), values = hue_pal()(12)[c(1, 7, 3, 9, 5)]) +
labs(title = "Nipbl Flox/+", x = "t-SNE 1", y = "t-SNE 2", color = "Embryo") +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = tsne.flox, filename = "tsne-by-embryo-flox-all-aggr.png", device = "png", width = 6.5, height = 4.5)

tsne.fin <- ggplot() +
geom_point(data = sample(tsne.whole %>% filter(genotype == "FIN")), mapping = aes(x = tsne.1, y = tsne.2, color = embryo.id), size = 0.1) +
scale_color_manual(labels = c("EP6", "EP7", "EP8", "EP9", "EP10", "EP11"), values = hue_pal()(12)[c(2, 8, 4, 10, 6, 12)]) +
labs(title = "Nipbl FIN/+", x = "t-SNE 1", y = "t-SNE 2", color = "Embryo") +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = tsne.fin, filename = "tsne-by-embryo-fin-all-aggr.png", device = "png", width = 6.5, height = 4.5)

all.aggr <- FindNeighbors(all.aggr, dims = 1:pc.number)
save(all.aggr, file = "all-aggr-neighbors.Robj")

dir.create("all-aggr-cluster-iters")
setwd("/dfs6/pub/schea2/20210505-seurat/all-aggr-cluster-iters")

res = 0.1

for (i in 2:15) {
    
all.aggr <- FindClusters(all.aggr, resolution = res, verbose = TRUE)

cluster.number = i
step = 0.1
loop.count = 1
loop.vector = vector()

while (max(as.numeric(all.aggr@active.ident[])) != cluster.number & (loop.count < 1000))
{
    loop.count = loop.count + 1
    
    #First, check to see if we're in stuck in a loop, and if we are, decrement the STEP
    if(length(loop.vector) > 3)
    {
        recent <- loop.vector[length(loop.vector)]
        second <- loop.vector[length(loop.vector) - 1]
        third <- loop.vector[length(loop.vector) - 2]
        if((recent == third) & (recent != second))
        {
            step = step/10
            print(paste0("step = ", step))
        }
    }
    #If the max cluster # is less than 15, we need a slightly higher resolution.
    if(max(as.numeric(all.aggr@active.ident[])) < cluster.number)
    {
        res = res + step
        all.aggr <- FindClusters(object = all.aggr, resolution = res, verbose = TRUE)
        loop.vector <- c(loop.vector, max(as.numeric(all.aggr@active.ident[])))
    }
    
    #If the max cluster is GREATER than 15, we need a slightly lower resolution
    if(max(as.numeric(all.aggr@active.ident[])) > cluster.number)
    {
        res = res - step
        all.aggr <- FindClusters(object = all.aggr, resolution = res, verbose = TRUE)
        loop.vector <- c(loop.vector, max(as.numeric(all.aggr@active.ident[])))
    }
    
    print(paste0("res = ", res))
}

res
length(levels(all.aggr))

renumbered.cluster.ids <- 1:i
names(renumbered.cluster.ids) <- levels(all.aggr)
all.aggr <- RenameIdents(all.aggr, renumbered.cluster.ids)
all.aggr.active.ident <- all.aggr@active.ident
filename.robj <- paste0("cluster-per-cell-all-aggr-clust-", i, ".Robj")
save(all.aggr.active.ident, file = filename.robj)

whole.tsne <- data.frame(all.aggr@reductions$tsne@cell.embeddings, all.aggr@active.ident)
colnames(whole.tsne) <- c("tsne.1", "tsne.2", "cluster")

whole.label <- whole.tsne %>% group_by(cluster) %>% summarize(tsne.1 = median(tsne.1), tsne.2 = median(tsne.2))

cluster.tsne <- ggplot() +
geom_point(data = sample(whole.tsne), mapping = aes(x = tsne.1, y = tsne.2, color = cluster), size = 0.1) +
geom_text(data = whole.label, mapping = aes(x = tsne.1, y = tsne.2, label = cluster), size = 3.5) +
labs(title = "Whole, Nipbl Flox/+ & FIN/+", x = "t-SNE 1", y = "t-SNE 2", color = "Cluster") +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

filename.png <- paste0("tsne-by-cluster-all-aggr-clust-", i, ".png")
ggsave(plot = cluster.tsne, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

#heatmap must be done on local computer
res = res + step
}

