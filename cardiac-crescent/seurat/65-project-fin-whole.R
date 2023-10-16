print("loading libraries")
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)

print("loading flox cells")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-neighbors.Robj")

print("loading flox pc's")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/sig-pc-flox-integr.Robj")

print("loading flox cluster identities")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/cluster-per-cell-flox-integr-res-0.0120.Robj")

print("setting flox cluster identities")
flox.integr$cluster.whole <- cluster.whole
Idents(flox.integr) <- cluster.whole

print("loading all cells")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/all-aggr-neighbors.Robj")
Idents(all.aggr) <- all.aggr$genotype

print("subsetting fin cells")
fin.aggr <- subset(all.aggr, idents = "Nipbl FIN/+")

remove(all.aggr)

print("performing sctransform on all fin cells")
fin.aggr <- SCTransform(object = fin.aggr, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.aggr$label.embryo))), n_genes = round(mean(fin.aggr$nFeature_RNA)), new.assay.name = "SCT_GENOTYPE", verbose = TRUE)

Idents(fin.aggr) <- fin.aggr$label.embryo

print("splitting fin cells by embryo")
embryo.fin <- SplitObject(fin.aggr, split.by = "label.embryo")

remove(fin.aggr)

print("performing sc transform on each fin embryo")
for (i in 1:length(embryo.fin)) {
  
  embryo.fin[[i]] <- SCTransform(embryo.fin[[i]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = TRUE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.fin[[i]]$nFeature_RNA)), verbose = TRUE)
  
  DefaultAssay(embryo.fin[[i]]) <- "SCT_EMBRYO"
  
}

print("grabbing variable genes from flox")
all.var.genes <- flox.integr@assays$integrated@var.features
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

print("loading fin cluster identities")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151/cluster-per-cell-fin-integr-res-0.0151.Robj")

load("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151/umap-by-cluster-fin-integr-res-0.0151-11.Robj")
cluster.clust <- umap.whole$cluster.whole

load("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151/umap-by-cluster-fin-integr-res-0.0151-28.Robj")
cluster.subclust <- umap.whole$cluster.whole

print("setting fin cluster identities")
fin.integr$cluster.subclust.fin <- cluster.subclust
fin.integr$cluster.clust.fin <- cluster.clust
fin.integr$cluster.whole.fin <- cluster.whole
Idents(fin.integr) <- fin.integr$cluster.whole.fin

print("finding transfer anchors")
fin.anchors <- FindTransferAnchors(reference = flox.integr, query = fin.integr, normalization.method = "SCT", reference.assay = "integrated", query.assay = "integrated",  reduction = "pcaproject", project.query = FALSE, reference.reduction = "pca", npcs = NULL, dims = 1:pc.number)
    
print("calculating projection scores")
fin.projected <- TransferData(anchorset = fin.anchors, refdata = flox.integr@active.ident, dims = 1:pc.number)

print("saving projection scores")
fin.projected$predicted.id <- factor(fin.projected$predicted.id, levels = levels(flox.integr$cluster.whole))

dir.create("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120")

save(fin.projected, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projection-scores-fin-integr-projected-flox-integr-res-0.0120.Robj")
  
print("integrating embeddings")
fin.integr <- IntegrateEmbeddings(anchorset = fin.anchors, new.reduction.name = "ref.pca", reference = flox.integr, query = fin.integr)

print("projecting umap")
fin.integr <- ProjectUMAP(query = fin.integr, query.dims = 1:pc.number, reference = flox.integr, reference.dims = 1:pc.number, query.reduction = "ref.pca", reduction.model = "umap", reference.reduction = "pca")

print("saving cluster identities")
cluster.whole.flox <- fin.projected$predicted.id
save(cluster.whole.flox, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-cluster-per-cell-fin-integr-projected-flox-integr-res-0.0120.Robj")

print("setting fin cluster identities")
fin.integr$cluster.whole.flox <- cluster.whole.flox
Idents(fin.integr) <- cluster.whole.flox

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

print("grabbing umap coordinates")
whole.umap <- data.frame(fin.integr@reductions$ref.umap@cell.embeddings, fin.integr$stage, fin.integr$label.embryo, fin.integr$cluster.whole.fin, fin.integr$cluster.clust.fin, fin.integr$cluster.subclust.fin, fin.integr$cluster.whole.flox)

colnames(whole.umap) <- c("umap.1", "umap.2", "stage", "embryo", "cluster.whole.fin", "cluster.clust.fin", "cluster.subclust.fin", "cluster.whole.flox")

print("saving umap coordinates")
save(whole.umap, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-by-stage-embryo-cluster-fin-integr-projected-flox-integr-res-0.0120.Robj")

stage.umap <- ggplot() +
  geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
  labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = stage.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-by-stage-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 3.5, height = 2.25)

stage.umap <- ggplot() +
  geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
  facet_wrap(~stage, ncol = 2) +
  labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

ggsave(plot = stage.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-facet-stage-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 6.5, height = 2.25)

embryo.umap <- ggplot() +
  geom_point(data = sample(whole.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
  scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
  labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = embryo.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-by-embryo-late-bud-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 3.5, height = 2.25)

embryo.umap <- ggplot() +
  geom_point(data = sample(whole.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
  facet_wrap(~embryo, ncol = 4) +
  scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
  labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

ggsave(plot = embryo.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-facet-embryo-late-bud-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 6.5, height = 3.6)

embryo.umap <- ggplot() +
  geom_point(data = sample(whole.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
  scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
  labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

ggsave(plot = embryo.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-by-embryo-cardiac-crescent-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 3.5, height = 2.25)

embryo.umap <- ggplot() +
  geom_point(data = sample(whole.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
  facet_wrap(~embryo, ncol = 4) +
  scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
  labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

ggsave(plot = embryo.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-facet-embryo-cardiac-crescent-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 6.5, height = 3.6)

whole.label <- whole.umap %>% group_by(cluster.whole.fin) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() + 
  geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.whole.fin), shape = 20, stroke = 0, size = 0.5) +
  geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster.whole.fin), size = 2.3) +
  labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Nipbl FIN/+\nCluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-by-fin-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 3.5, height = 2.25, units = "in")

whole.label <- whole.umap %>% group_by(cluster.whole.flox) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

cluster.umap <- ggplot() + 
  geom_point(data = sample(whole.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.whole.flox), shape = 20, stroke = 0, size = 0.5) +
  geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster.whole.flox), size = 2.3) +
  labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Nipbl Flox/+\nCluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = cluster.umap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-umap-by-flox-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", width = 3.5, height = 2.25, units = "in")

rm(flox.integr)
rm(fin.anchors)
rm(fin.projected)

print("saving fin projected")
save(fin.integr, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/-fin-integr-projected-flox-integr-res-0.0120.Robj")

print("preparing to find DEGs")
fin.integr <- PrepSCTFindMarkers(fin.integr, assay = "SCT_GENOTYPE", verbose = TRUE)

print("finding DEGs")
diff.genes <- FindAllMarkers(fin.integr, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.integr)[2]/length(levels(fin.integr$cluster.whole.flox)), return.thresh = 0.05)

print("saving DEGs")
save(diff.genes, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-diff-genes-fin-integr-projected-flox-integr-res-0.0120.Robj")

write.csv(x = diff.genes, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/projected-diff-genes-fin-integr-projected-flox-integr-res-0.0120.csv")

top.genes <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(fin.integr$cluster.whole.flox)))

cell.names <- rownames(fin.integr@meta.data)

heatmap <- DoHeatmap(fin.integr, cells = sample(cell.names, round(length(cell.names)*0.5), replace = FALSE), assay = "SCT_GENOTYPE", slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
  scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
  scale_x_discrete(position = "top") +
  labs(title = "Nipbl FIN/+", x = "Nipbl Flox/+ Cluster", y = "Top Genes", fill = "Standardized\nLog1p UMI") +
  guides(color = FALSE) +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

ggsave(plot = heatmap, filename = "/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/heatmap-diff-genes-fin-integr-projected-flox-integr-res-0.0120.png", device = "png", units = "in", width = 6.5, height = 9)
