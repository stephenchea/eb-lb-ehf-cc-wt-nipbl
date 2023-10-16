library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(pracma)
library(patchwork)

setwd("/dfs6/pub/schea2/20210505-seurat")

load("/dfs6/pub/schea2/20210505-seurat/all-aggr-neighbors.Robj")
all.aggr$genotype <- factor(all.aggr$genotype, levels = c("Flox", "FIN"))
Idents(all.aggr) <- all.aggr$genotype

#subset flox cells from all cells
flox.aggr <- subset(all.aggr, idents = "Flox")

flox.aggr <- SCTransform(object = flox.aggr, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = round(dim(flox.aggr)[2]*.5), n_genes = round(dim(flox.aggr)[1]*.1), new.assay.name = "SCT_GENOTYPE", verbose = TRUE)

#split flox cells according to sample id
flox.embryo <- SplitObject(flox.aggr, split.by = "embryo.id")

for (i in 1:length(flox.embryo)) {
    
    flox.embryo[[i]] <- SCTransform(flox.embryo[[i]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = NULL, n_genes = round(dim(flox.embryo[[i]])[1]*.1), verbose = TRUE)
    
    DefaultAssay(flox.embryo[[i]]) <- "SCT_EMBRYO"
}

print("grabbing variable genes from each embryo")
var.genes <- list(flox.embryo[[1]]@assays$SCT_EMBRYO@var.features,
flox.embryo[[2]]@assays$SCT_EMBRYO@var.features,
flox.embryo[[3]]@assays$SCT_EMBRYO@var.features,
flox.embryo[[4]]@assays$SCT_EMBRYO@var.features,
flox.embryo[[5]]@assays$SCT_EMBRYO@var.features
)

names(var.genes) <- names(flox.embryo)

lapply(var.genes, FUN = length)

print("finding only those variable genes that are unique")
uniq.genes <- unique(c(var.genes[[1]],
var.genes[[2]],
var.genes[[3]],
var.genes[[4]],
var.genes[[5]]
))

print("counting how many times each gene occurs")
count.gene <- list()
for (e in uniq.genes) {
    for(r in names(var.genes)) {
        count.gene[[e]][[r]]<- e %in% var.genes[[r]]
    }
}

print("ranking unique variable genes")
count.rank <- lapply(X = count.gene, FUN = unlist)
count.rank <- lapply(X = count.rank, FUN = sum)
count.rank <- sort(unlist(count.rank), decreasing = TRUE)

print("table of rankings")
print(table(count.rank))
print("finding minimum")
print(which.min(table(count.rank)))
print("how many genes")
print(sum(count.rank >= which.min(table(count.rank))))

print("storing number of genes for integration")
no.genes <- sum(count.rank >= which.min(table(count.rank)))

print("grabbing genes for integration")
feat.integr <- SelectIntegrationFeatures(object.list = flox.embryo, nfeatures = no.genes)

flox.embryo <- PrepSCTIntegration(object.list = flox.embryo, anchor.features = feat.integr)

flox.anchors <- FindIntegrationAnchors(object.list = flox.embryo, normalization.method = "SCT", anchor.features = feat.integr)

#find genes common to all Flox embryos
com.genes <- Reduce(intersect, list(rownames(flox.embryo[[1]]@assays$SCT_EMBRYO@data),
rownames(flox.embryo[[2]]@assays$SCT_EMBRYO@data),
rownames(flox.embryo[[3]]@assays$SCT_EMBRYO@data),
rownames(flox.embryo[[4]]@assays$SCT_EMBRYO@data),
rownames(flox.embryo[[5]]@assays$SCT_EMBRYO@data)
))

length(com.genes)

flox.integr <- IntegrateData(anchorset = flox.anchors, normalization.method = "SCT", features.to.integrate = com.genes)

DefaultAssay(flox.integr) <- "integrated"

print("performing principal component reduction")
flox.integr <- RunPCA(flox.integr, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

print("grabbing deviations of principal components")
stdev.pca <- flox.integr@reductions$pca@stdev

var.exp <- stdev.pca^2/sum(stdev.pca^2)

print("set number of principal components")
d.var.exp <- var.exp - median(var.exp)
pc.number <- sum(d.var.exp > mad(var.exp)*2)
save(pc.number, file = "sig-pc-flox-integr.Robj")

flox.integr <- RunUMAP(flox.integr, reduction = "pca", dims = 1:pc.number)

whole.umap <- data.frame(flox.integr@reductions$umap@cell.embeddings, flox.integr$embryo.id)
colnames(whole.umap) <- c("umap.1", "umap.2", "embryo.id")
save(whole.umap, file = "umap-by-embryo-flox-integr-umap.Robj")

umap.embryo <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = embryo.id), size = 0.1) +
scale_color_manual(labels = c("EP1", "EP2", "EP3", "EP4", "EP5"), values = hue_pal()(12)[c(1, 7, 3, 9, 5)]) +
labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo\nPair") +
guides(color = guide_legend(override.aes = list(size = 3.5))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = umap.embryo, filename = "umap-by-embryo-flox-integr.png", device = "png", width = 6.5, height = 4.5)

flox.integr <- FindNeighbors(flox.integr, dims = 1:pc.number)
save(flox.integr, file = "flox-integr-neighbors.Robj")

dir.create("flox-integr-cluster-iters")
setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-cluster-iters")

res = 0.1

for (i in 2:15) {
    
flox.integr <- FindClusters(flox.integr, resolution = res, verbose = TRUE)

cluster.number = i
step = 0.1
loop.count = 1
loop.vector = vector()

while (max(as.numeric(flox.integr@active.ident[])) != cluster.number & (loop.count < 1000))
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
    if(max(as.numeric(flox.integr@active.ident[])) < cluster.number)
    {
        res = res + step
        flox.integr <- FindClusters(object = flox.integr, resolution = res, verbose = TRUE)
        loop.vector <- c(loop.vector, max(as.numeric(flox.integr@active.ident[])))
    }
    
    #If the max cluster is GREATER than 15, we need a slightly lower resolution
    if(max(as.numeric(flox.integr@active.ident[])) > cluster.number)
    {
        res = res - step
        flox.integr <- FindClusters(object = flox.integr, resolution = res, verbose = TRUE)
        loop.vector <- c(loop.vector, max(as.numeric(flox.integr@active.ident[])))
    }
    
    print(paste0("res = ", res))
}

res
length(levels(flox.integr))

renumbered.cluster.ids <- 1:i
names(renumbered.cluster.ids) <- levels(flox.integr)
flox.integr <- RenameIdents(flox.integr, renumbered.cluster.ids)
flox.integr.active.ident <- flox.integr@active.ident
filename.robj <- paste0("cluster-per-cell-flox-integr-clust-", i, ".Robj")
save(flox.integr.active.ident, file = filename.robj)

whole.umap <- data.frame(flox.integr@reductions$umap@cell.embeddings, flox.integr@active.ident)
colnames(whole.umap) <- c("umap.1", "umap.2", "cluster")

umap.label <- whole.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 3.5) +
labs(title = "Whole, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 3.5))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-cluster-flox-integr-clust-", i, ".png")
ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

#heatmap must be done on local computer
res = res + step
}
