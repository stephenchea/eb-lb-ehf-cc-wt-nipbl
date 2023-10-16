library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(pracma)

setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

load("/share/crsp/lab/alcalof/schea2/20220414-seurat/all-aggr-neighbors.Robj")
Idents(all.aggr) <- all.aggr$genotype

#subset fin cells from all cells
fin.aggr <- subset(all.aggr, idents = "Nipbl FIN/+")

remove(all.aggr)

fin.aggr <- SCTransform(object = fin.aggr, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.aggr$label.embryo))), n_genes = round(mean(fin.aggr$nFeature_RNA)), new.assay.name = "SCT_GENOTYPE", verbose = TRUE)

Idents(fin.aggr) <- fin.aggr$label.embryo

#split fin cells according to sample id
embryo.fin <- SplitObject(fin.aggr, split.by = "label.embryo")

remove(fin.aggr)

for (i in 1:length(embryo.fin)) {
    
    embryo.fin[[i]] <- SCTransform(embryo.fin[[i]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = TRUE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.fin[[i]]$nFeature_RNA)), verbose = TRUE)
    
    DefaultAssay(embryo.fin[[i]]) <- "SCT_EMBRYO"
    
    }

print("grabbing variable genes from each embryo")
var.genes <- list(embryo.fin[[1]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[2]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[3]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[4]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[5]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[6]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[7]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[8]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[9]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[10]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[11]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[12]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[13]]@assays$SCT_EMBRYO@var.features,
embryo.fin[[14]]@assays$SCT_EMBRYO@var.features)

names(var.genes) <- names(embryo.fin)

print("number of variable genes per embryo")
lapply(var.genes, FUN = length)

print("finding unique variable genes among late bud embryos")
lb.uniq.genes <- unique(c(var.genes[[1]],
var.genes[[2]],
var.genes[[3]],
var.genes[[4]],
var.genes[[5]],
var.genes[[6]]
))

print("number of unique variable genes in late bud embyros")
length(lb.uniq.genes)

print("counting how many times each gene occurs")
lb.count.gene <- list()
for (e in lb.uniq.genes) {
    for(r in names(var.genes)[1:6]) {
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
lb.var.genes <- names(lb.count.rank[lb.count.rank >= 6])

print("finding unique genes among cardiac crescent embryos")
cc.uniq.genes <- unique(c(var.genes[[7]],
var.genes[[8]],
var.genes[[9]],
var.genes[[10]],
var.genes[[11]],
var.genes[[12]],
var.genes[[13]],
var.genes[[14]]
))

print("number of unique variable genes in cardiac crescent embyros")
length(cc.uniq.genes)

print("counting how many times each gene occurs")
cc.count.gene <- list()
for (e in cc.uniq.genes) {
    for(r in names(var.genes)[7:14]) {
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
sct.genes <- list(rownames(embryo.fin[[1]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[2]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[3]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[4]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[5]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[6]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[7]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[8]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[9]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[10]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[11]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[12]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[13]]@assays$SCT_EMBRYO),
rownames(embryo.fin[[14]]@assays$SCT_EMBRYO))

print("keep only variable genes among common sct genes")
feat.integr <- intersect(all.var.genes, Reduce(intersect, sct.genes))
length(feat.integr)

embryo.fin <- PrepSCTIntegration(object.list = embryo.fin, anchor.features = feat.integr)

for (i in 1:length(embryo.fin)) {
    embryo.fin[[i]] <- RunPCA(embryo.fin[[i]], assay = "SCT_EMBRYO", features = feat.integr, npcs = 50, verbose = TRUE)
}

anchors.fin <- FindIntegrationAnchors(object.list = embryo.fin, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca")

fin.integr <- IntegrateData(anchorset = anchors.fin, normalization.method = "SCT", features.to.integrate = feat.integr)

DefaultAssay(fin.integr) <- "integrated"

fin.integr[["SCT_EMBRYO"]] <- NULL

print("performing principal component reduction")
fin.integr <- RunPCA(fin.integr, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

print("grabbing deviations of principal components")
stdev.pca <- fin.integr@reductions$pca@stdev
dens.stdev <- density(stdev.pca)
peak.stdev <- findpeaks(x = dens.stdev$y*-1, npeaks = 1)[1, 2]
elbow.stdev <- dens.stdev$x[peak.stdev]
pc.number <- sum(stdev.pca > elbow.stdev)
print(pc.number)

save(pc.number, file = "sig-pc-fin-integr.Robj")

fin.integr <- RunUMAP(fin.integr, reduction = "pca", dims = 1:pc.number, return.model = TRUE)

fin.integr$stage <- factor(fin.integr$stage, levels = c("LB", "CC"))
fin.integr$label.embryo <- factor(fin.integr$label.embryo, levels = c("LB6",
"LB7",
"LB8",
"LB9",
"LB10",
"LB11",
"CC9",
"CC10",
"CC11",
"CC12",
"CC13",
"CC14",
"CC15",
"CC16"))

whole.umap <- data.frame(fin.integr@reductions$umap@cell.embeddings, fin.integr$stage, fin.integr$label.embryo)
colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "embryo")
save(whole.umap, file = "umap-by-embryo-fin-integr-umap.Robj")

stage.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
    labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

#ggsave(plot = stage.umap, filename = "umap-by-stage-fin-integr.png", device = "png", width = 3.5, height = 2.25)

stage.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~stage, ncol = 2) +
    labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

#ggsave(plot = stage.umap, filename = "umap-facet-stage-fin-integr.png", device = "png", width = 6.5, height = 2.25)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
    labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

#ggsave(plot = embryo.umap, filename = "umap-by-embryo-late-bud-fin-integr.png", device = "png", width = 3.5, height = 2.25)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~embryo, ncol = 4) +
    scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
    labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

#ggsave(plot = embryo.umap, filename = "umap-facet-embryo-late-bud-fin-integr.png", device = "png", width = 6.5, height = 3.6)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
    labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

#ggsave(plot = embryo.umap, filename = "umap-by-embryo-cardiac-crescent-fin-integr.png", device = "png", width = 3.5, height = 2.25)

embryo.umap <- ggplot() +
geom_point(data = sample(whole.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~embryo, ncol = 4) +
    scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
    labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

#ggsave(plot = embryo.umap, filename = "umap-facet-embryo-cardiac-crescent-fin-integr.png", device = "png", width = 6.5, height = 3.6)

print("finding neighbors among fin cells")
fin.integr <- FindNeighbors(fin.integr, dims = 1:pc.number)
save(fin.integr, file = "fin-integr-neighbors.Robj")

dir.create("fin-integr-cluster-iters")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-cluster-iters")

for (i in seq(from = 0, to = 0.5, by = 0.0001)) {
    
fin.integr <- FindClusters(fin.integr, resolution = i, verbose = TRUE)

renumbered.cluster.ids <- 1:length(levels(fin.integr))
names(renumbered.cluster.ids) <- levels(fin.integr)
fin.integr <- RenameIdents(fin.integr, renumbered.cluster.ids)
fin.integr.active.ident <- fin.integr@active.ident
filename.robj <- paste0("cluster-per-cell-fin-integr-res-", sprintf(fmt = "%.4f", i), ".Robj")
save(fin.integr.active.ident, file = filename.robj)

whole.umap <- data.frame(fin.integr@reductions$umap@cell.embeddings, fin.integr$stage, fin.integr$label.embryo, fin.integr@active.ident)
colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "embryo", "cluster")

whole.label <- whole.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() +
geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), shape = 20, stroke = 0, size = 0.5) +
geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 2.3) +
labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
guides(color = guide_legend(override.aes = list(size = 2.3))) +
theme_classic(base_size = 7) +
theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

filename.png <- paste0("umap-by-cluster-fin-integr-res-", sprintf(fmt = "%.4f", i), ".png")
##ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

}
