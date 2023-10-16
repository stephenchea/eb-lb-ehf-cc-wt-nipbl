library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

#set working directory
setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18")

#load seurat object
load("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18/fin-integr-projected-flox-integr-clust-5-18.Robj")

#load wildtype identities
load("projected-umap-by-wildtype-fin-integr-projected-flox-integr-clust-5-18.Robj")
fin.integr$wildtype.whole <- whole.umap$wildtype.whole

#load germ layer identities
load("projected-umap-by-germ-layer-fin-integr-projected-flox-integr-clust-5-18.Robj")
fin.integr$germ.whole <- whole.umap$germ.whole

#set identity to germ layer
Idents(fin.integr) <- whole.umap$germ.whole

#begin loop across germ layers
for (i in c("ECT")) {
    
print(i)
print("subsetting cluster")
fin.germ <- subset(fin.integr, idents = i)

print("setting identities to cluster")
fin.germ$wildtype.whole <- factor(fin.germ$wildtype.whole)
Idents(fin.germ) <- fin.germ$wildtype.whole

print("creating new directory")
direct.name <- paste0("fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)))
dir.create(direct.name)

print("set new directory")
direct.path <- paste0("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18/", direct.name)
setwd(direct.path)

print("applying sc transform to entire cluster")
assay.name <- paste0("SCT_", i)
fin.germ <- SCTransform(object = fin.germ, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = round(dim(fin.germ)[2]*.5), n_genes = round(dim(fin.germ)[1]*.1), new.assay.name = assay.name, verbose = TRUE)

print("splitting by embryo")
fin.embryo <- SplitObject(fin.germ, split.by = "embryo.id")

#open bracket 2
print("applying sc transform across embryos")
for (t in 1:length(fin.embryo)) {
    fin.embryo[[t]] <- SCTransform(fin.embryo[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = NULL, n_genes = round(dim(fin.embryo[[t]])[1]*.1), verbose = TRUE)
    DefaultAssay(fin.embryo[[t]]) <- "SCT_EMBRYO"
}
#close open bracket 2

print("grabbing variable genes from each embryo")
var.genes <- list(fin.embryo[[1]]@assays$SCT_EMBRYO@var.features,
fin.embryo[[2]]@assays$SCT_EMBRYO@var.features,
fin.embryo[[3]]@assays$SCT_EMBRYO@var.features,
fin.embryo[[4]]@assays$SCT_EMBRYO@var.features,
fin.embryo[[5]]@assays$SCT_EMBRYO@var.features,
fin.embryo[[6]]@assays$SCT_EMBRYO@var.features
)
names(var.genes) <- names(fin.embryo)
lapply(var.genes, FUN = length)

print("finding only those variable genes that are unique")
uniq.genes <- unique(c(var.genes[[1]],
var.genes[[2]],
var.genes[[3]],
var.genes[[4]],
var.genes[[5]],
var.genes[[6]]
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
feat.integr <- SelectIntegrationFeatures(object.list = fin.embryo, nfeatures = no.genes)

print("performing principal component reduction")
fin.germ <- RunPCA(fin.germ, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

print("grabbing deviations of principal components")
stdev.pca <- fin.germ@reductions$pca@stdev
var.exp <- stdev.pca^2/sum(stdev.pca^2)

print("set number of principal components")
d.var.exp <- var.exp - median(var.exp)
pc.number <- sum(d.var.exp > mad(var.exp)*2)

filename.robj <- paste0("sig-pc-fin-", i, ".Robj")
save(object = pc.number, file = filename.robj)

fin.germ <- RunUMAP(fin.germ, reduction = "pca", dims = 1:pc.number, return.model = TRUE)

filename.robj <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), "/flox-", i, "-umap.Robj")
load(file = filename.robj)

filename.robj <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), "/sig-pc-flox-", i, ".Robj")
load(file = filename.robj)

fin.anchors <- FindTransferAnchors(reference = flox.germ, query = fin.germ, normalization.method = "SCT", reference.assay = "integrated", query.assay = "integrated",  reduction = "pcaproject", project.query = FALSE, reference.reduction = "pca", npcs = NULL, dims = 1:pc.number, max.features = length(rownames(flox.germ))*0.1)

fin.germ <- IntegrateEmbeddings(anchorset = fin.anchors, new.reduction.name = "ref.pca", reference = flox.germ, query = fin.germ)

fin.germ <- ProjectUMAP(query = fin.germ, query.dims = 1:pc.number, reference = flox.germ, reference.dims = 1:pc.number,
    query.reduction = "ref.pca", reference.reduction = "pca", reduction.model = "umap")
    
filename.robj <- paste0("fin-", i, "-umap.Robj")
save(fin.germ, file = filename.robj)
#########
germ.umap <- data.frame(fin.germ@reductions$ref.umap@cell.embeddings, fin.germ$embryo.id)
colnames(germ.umap) <- c("umap.1", "umap.2", "embryo.id")

filename.robj <- paste0("umap-by_embryo-fin-", i, "-umap.Robj")
save(germ.umap, file = filename.robj)

umap.embryo <- ggplot() +
geom_point(data = sample(germ.umap), mapping = aes(x = umap.1, y = umap.2, color = embryo.id), size = 0.1) +
scale_color_manual(labels = c("EP6", "EP7", "EP8", "EP9", "EP10", "EP11"), values = hue_pal()(12)[c(2, 8, 4, 10, 6, 12)]) +
labs(title = paste0(i, ", Nipbl FIN/+"), x = "UMAP 1", y = "UMAP 2", color = "Embryo\nPair") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-embryo-fin-", i, ".png")
ggsave(plot = umap.embryo, filename = filename.png, device = "png", width = 6.5, height = 4.5)

fin.germ.active.ident <- fin.germ$wildtype.whole
filename.robj <- paste0("cluster-per-cell-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".Robj")
save(fin.germ.active.ident, file = filename.robj)

germ.umap <- data.frame(fin.germ@reductions$ref.umap@cell.embeddings, fin.germ$embryo.id, fin.germ$wildtype.whole)
colnames(germ.umap) <- c("umap.1", "umap.2", "embryo.id", "wildtype.whole")

filename.robj <- paste0("umap-by-cluster-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".Robj")
save(germ.umap, file = filename.robj)

umap.label <- germ.umap %>% group_by(wildtype.whole) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(germ.umap), mapping = aes(x = umap.1, y = umap.2, color = wildtype.whole), size = 0.1) +
geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = wildtype.whole), size = 3.5) +
labs(title = paste0(i, ", Nipbl FIN/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-cluster-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".png")
ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

marker.genes <- FindAllMarkers(fin.germ, assay = assay.name, slot = "data", only.pos = TRUE, test.use = "roc", min.pct = 0.2, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, verbose = TRUE)

filename.robj <- paste0("marker-genes-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".Robj")
save(marker.genes, file = filename.robj)

filename.csv <- paste0("marker-genes-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".csv")
write.csv(x = marker.genes, file = filename.csv)

top.markers <- marker.genes %>% group_by(cluster) %>% slice_max(myAUC, n = 71/length(levels(fin.germ$wildtype.whole)))

heatmap <- DoHeatmap(fin.germ, assay = assay.name, slot = "scale.data", features = as.vector(top.markers$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = paste0(i, ", Nipbl FIN/+"), x = "Cluster", y = "Top Genes by AUC", fill = "Centered\nNorm UMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("heatmap-marker-genes-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18")

}
