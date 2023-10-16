library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(pracma)

print("setting working directory")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

print("load flox cells after sctranform")
load(file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-aggr-sctrans.Robj")

print("grabbing variable genes from each embryo")
var.genes <- list(embryo.flox[[1]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[2]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[3]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[4]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[5]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[6]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[7]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[8]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[9]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[10]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[11]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[12]]@assays$SCT_EMBRYO@var.features,
embryo.flox[[13]]@assays$SCT_EMBRYO@var.features)

names(var.genes) <- names(embryo.flox)

print("number of variable genes per embryo")
lapply(var.genes, FUN = length)

print("finding unique variable genes among late bud embryos")
lb.uniq.genes <- unique(c(var.genes[[1]],
var.genes[[2]],
var.genes[[3]],
var.genes[[4]],
var.genes[[5]]
))

print("number of unique variable genes in late bud embyros")
length(lb.uniq.genes)

print("counting how many times each gene occurs")
lb.count.gene <- list()
for (e in lb.uniq.genes) {
    for(r in names(var.genes)[1:5]) {
        lb.count.gene[[e]][[r]]<- e %in% var.genes[[r]]
    }
}

print("ranking unique variable genes")
lb.count.rank <- lapply(X = lb.count.gene, FUN = unlist)
lb.count.rank <- lapply(X = lb.count.rank, FUN = sum)
lb.count.rank <- sort(unlist(lb.count.rank), decreasing = TRUE)

print("table of rankings")
print(table(lb.count.rank))

print("grabbing variable genes common to all late bud embryos")
lb.var.genes <- names(lb.count.rank[lb.count.rank >= 5])

print("finding unique genes among cardiac crescent embryos")
cc.uniq.genes <- unique(c(var.genes[[6]],
var.genes[[7]],
var.genes[[8]],
var.genes[[9]],
var.genes[[10]],
var.genes[[11]],
var.genes[[12]],
var.genes[[13]]
))

print("number of unique variable genes in cardiac crescent embyros")
length(cc.uniq.genes)

print("counting how many times each gene occurs")
cc.count.gene <- list()
for (e in cc.uniq.genes) {
    for(r in names(var.genes)[6:13]) {
        cc.count.gene[[e]][[r]]<- e %in% var.genes[[r]]
    }
}

print("ranking unique variable genes")
cc.count.rank <- lapply(X = cc.count.gene, FUN = unlist)
cc.count.rank <- lapply(X = cc.count.rank, FUN = sum)
cc.count.rank <- sort(unlist(cc.count.rank), decreasing = TRUE)

print("table of rankings")
print(table(cc.count.rank))

print("grabbing variable genes common to all cardiac crescent embryos")
cc.var.genes <- names(cc.count.rank[cc.count.rank >= 8])

print("grabbing variable genes common to late bud and cardiac crescent embryos")
all.var.genes <- unique(c(lb.var.genes, cc.var.genes))

print("removing mitochondrial genes from variable genes")
all.var.genes <- str_subset(all.var.genes, "mt-", negate = TRUE)

length(all.var.genes)

print("grabbing sct genes from each embryo")
sct.genes <- list(rownames(embryo.flox[[1]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[2]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[3]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[4]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[5]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[6]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[7]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[8]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[9]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[10]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[11]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[12]]@assays$SCT_EMBRYO),
rownames(embryo.flox[[13]]@assays$SCT_EMBRYO))

print("keep only variable genes among common sct genes")
feat.integr <- intersect(all.var.genes, Reduce(intersect, sct.genes))
length(feat.integr)

embryo.flox <- PrepSCTIntegration(object.list = embryo.flox, anchor.features = feat.integr)

for (i in 1:length(embryo.flox)) {
    embryo.flox[[i]] <- RunPCA(embryo.flox[[i]], assay = "SCT_EMBRYO", features = feat.integr, npcs = 50, verbose = TRUE)
}

anchors.flox <- FindIntegrationAnchors(object.list = embryo.flox, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca")

flox.integr <- IntegrateData(anchorset = anchors.flox, normalization.method = "SCT", features.to.integrate = feat.integr)

DefaultAssay(flox.integr) <- "integrated"

flox.integr[["SCT_EMBRYO"]] <- NULL

print("performing principal component reduction")
flox.integr <- RunPCA(flox.integr, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

print("grabbing deviations of principal components")
stdev.pca <- flox.integr@reductions$pca@stdev
dens.stdev <- density(stdev.pca)
peak.stdev <- findpeaks(x = dens.stdev$y*-1, npeaks = 1)[1, 2]
elbow.stdev <- dens.stdev$x[peak.stdev]
pc.number <- sum(stdev.pca > elbow.stdev)
print(pc.number)

save(pc.number, file = "sig-pc-flox-integr.Robj")

flox.integr <- RunUMAP(flox.integr, reduction = "pca", dims = 1:pc.number, return.model = TRUE)

flox.integr$stage <- factor(flox.integr$stage, levels = c("LB", "CC"))
flox.integr$label.embryo <- factor(flox.integr$label.embryo, levels = c("LB1",
"LB2",
"LB3",
"LB4",
"LB5",
"CC1",
"CC2",
"CC3",
"CC4",
"CC5",
"CC6",
"CC7",
"CC8"))

whole.umap <- data.frame(flox.integr@reductions$umap@cell.embeddings, flox.integr$stage, flox.integr$label.embryo)
colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "embryo")
save(whole.umap, file = "umap-by-embryo-flox-integr-umap.Robj")

stage.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
    labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

#ggsave(plot = stage.umap, filename = "umap-by-stage-flox-integr.png", device = "png", width = 3.5, height = 2.25)

stage.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~stage, ncol = 2) +
    labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

#ggsave(plot = stage.umap, filename = "umap-facet-stage-flox-integr.png", device = "png", width = 6.5, height = 2.25)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    scale_color_manual(values = hue_pal()(28)[c(1, 15, 3, 17, 5)]) +
    labs(title = "Late Bud, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

#ggsave(plot = embryo.umap, filename = "umap-by-embryo-late-bud-flox-integr.png", device = "png", width = 3.5, height = 2.25)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~embryo, ncol = 4) +
    scale_color_manual(values = hue_pal()(28)[c(1, 15, 3, 17, 5)]) +
    labs(title = "Late Bud, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

#ggsave(plot = embryo.umap, filename = "umap-facet-embryo-late-bud-flox-integr.png", device = "png", width = 6.5, height = 3.6)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    scale_color_manual(values = hue_pal()(28)[c(25, 13, 27, 2, 16, 4, 18, 6)]) +
    labs(title = "Cardiac Crescent, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

#ggsave(plot = embryo.umap, filename = "umap-by-embryo-cardiac-crescent-flox-integr.png", device = "png", width = 3.5, height = 2.25)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~embryo, ncol = 4) +
    scale_color_manual(values = hue_pal()(28)[c(25, 13, 27, 2, 16, 4, 18, 6)]) +
    labs(title = "Cardiac Crescent, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

#ggsave(plot = embryo.umap, filename = "umap-facet-embryo-cardiac-crescent-flox-integr.png", device = "png", width = 6.5, height = 3.6)

print("finding neighbors among flox cells")
flox.integr <- FindNeighbors(flox.integr, dims = 1:pc.number)
save(flox.integr, file = "flox-integr-neighbors.Robj")

dir.create("flox-integr-cluster-iters")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-cluster-iters")

for (i in seq(from = 0, to = 0.5, by = 0.0001)) {
    
flox.integr <- FindClusters(flox.integr, resolution = i, verbose = TRUE)

renumbered.cluster.ids <- 1:length(levels(flox.integr))
names(renumbered.cluster.ids) <- levels(flox.integr)
flox.integr <- RenameIdents(flox.integr, renumbered.cluster.ids)
flox.integr.active.ident <- flox.integr@active.ident
filename.robj <- paste0("cluster-per-cell-flox-integr-res-", sprintf(fmt = "%.4f", i), ".Robj")
save(flox.integr.active.ident, file = filename.robj)

whole.umap <- data.frame(flox.integr@reductions$umap@cell.embeddings, flox.integr$stage, flox.integr$label.embryo, flox.integr@active.ident)
colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "embryo", "cluster")

whole.label <- whole.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), shape = 20, stroke = 0, size = 0.5) +
geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 2.3) +
labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-cluster-flox-integr-res-", sprintf(fmt = "%.4f", i), ".png")
##ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

}
