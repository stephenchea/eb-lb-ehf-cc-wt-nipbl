library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

res <- c("0.177")
names(res) <- c("1.1")

for (i in c("1.1")) {

    print("loading flox cells")
    filename.robj <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-1-res-0.074/flox-", i, "-neighbors.Robj")
    load(file = filename.robj)

    print("loading flox pc's")
    filename.robj <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-1-res-0.074/sig-pc-flox-", i, ".Robj")
    load(file = filename.robj)

    print("loading flox cluster identities")
    filename.robj <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-1-res-0.074/flox-", i, "-res-", res[i], "/cluster-per-cell-flox-", i, "-res-", res[i], ".Robj")
    load(file = filename.robj)

    print("setting flox cluster identities")
    flox.subclust$cluster.subclust <- cluster.subclust
    Idents(flox.subclust) <- cluster.subclust
    
print("loading whole fin projected flox whole")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/fin-1-projected-flox-1-res-0.074")
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/fin-1-projected-flox-1-res-0.074/fin-1-projected-flox-1-res-0.074.Robj")

DefaultAssay(fin.clust) <- "RNA"
fin.clust[["SCT_1"]] <- NULL
fin.clust[["INT_1"]] <- NULL

print(paste0("subsetting cluster ", i))
fin.subclust <- subset(fin.clust, idents = i)
remove(fin.clust)

print("generating assay names")
n <- str_replace_all(i,"[.]","_")
sct.assay.name <- paste0("SCT_", n)
int.assay.name <- paste0("INT_", n)

print("performing sctransform on entire cluster")
fin.subclust <- SCTransform(object = fin.subclust, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.subclust$label.embryo))), n_genes = round(mean(fin.subclust$nFeature_RNA)), new.assay.name = sct.assay.name, verbose = TRUE)

print("splitting by embryo")
embryo.fin <- SplitObject(fin.subclust, split.by = "label.embryo")

print("performing sc transform across embryos")
for (t in 1:length(embryo.fin)) {
    
    embryo.fin[[t]] <- SCTransform(embryo.fin[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = TRUE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.fin[[t]]$nFeature_RNA)), verbose = TRUE)
    
    DefaultAssay(embryo.fin[[t]]) <- "SCT_EMBRYO"
    
}

print("grabbing variable genes from flox")
all.var.genes <- flox.subclust@assays[[int.assay.name]]@var.features
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

print("keeping only variable genes among common sct genes")
feat.integr <- intersect(all.var.genes, Reduce(intersect, sct.genes))
length(feat.integr)

print("preparing sctransform assay for integration")
embryo.fin <- PrepSCTIntegration(object.list = embryo.fin, anchor.features = feat.integr)

min.cells <- min(table(fin.subclust$label.embryo))

print("performing pca across all embryos")
for (a in 1:length(embryo.fin)) {
    embryo.fin[[a]] <- RunPCA(embryo.fin[[a]], assay = "SCT_EMBRYO", features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)
}

print("finding integration anchors across all embryos")
anchors.fin <- FindIntegrationAnchors(object.list = embryo.fin, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca", k.score = min(c(30, min.cells*0.95)), k.filter = min(c(200, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

print("integrating embryos")
fin.subclust <- IntegrateData(anchorset = anchors.fin, normalization.method = "SCT", features.to.integrate = feat.integr, new.assay.name = int.assay.name, k.weight = min(c(100, min.cells*0.82)), dims = 1:min(c(30, min.cells*0.95)))

print("setting default assay to integrated assay")
DefaultAssay(fin.subclust) <- int.assay.name

fin.subclust[["SCT_EMBRYO"]] <- NULL

    print("finding transfer anchors")
    fin.anchors <- FindTransferAnchors(reference = flox.subclust, query = fin.subclust, normalization.method = "SCT", reference.assay = int.assay.name, query.assay = int.assay.name,  reduction = "pcaproject", project.query = FALSE, reference.reduction = "pca", npcs = NULL, dims = 1:pc.number)
    
    print("calculating projection scores")
    fin.projected <- TransferData(anchorset = fin.anchors, refdata = flox.subclust@active.ident, dims = 1:pc.number)

print("saving projection scores")
fin.projected$predicted.id <- factor(fin.projected$predicted.id, levels = levels(flox.subclust$cluster.subclust))

dir.name <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/fin-1-projected-flox-1-res-0.074/fin-", i, "-projected-flox-", i, "-res-", res[i])
dir.create(dir.name)

filename.robj <- paste0(dir.name, "/projection-scores-fin-", i, "-projected-flox-", i, "-res-", res[i], ".Robj")
save(fin.projected, file = filename.robj)

  print("integrating embeddings")
  fin.subclust <- IntegrateEmbeddings(anchorset = fin.anchors, new.reduction.name = "ref.pca", reference = flox.subclust, query = fin.subclust)

print("projecting umap")
fin.subclust <- ProjectUMAP(query = fin.subclust, query.dims = 1:pc.number, reference = flox.subclust, reference.dims = 1:pc.number, query.reduction = "ref.pca", reference.reduction = "pca", reduction.model = "umap")

print("saving cluster identities")
cluster.subclust.flox <- fin.projected$predicted.id
filename.robj <- paste0(dir.name, "/projected-cluster-per-cell-fin-", i, "-projected-flox-", i, "-res-", res[i], ".Robj")
save(cluster.subclust.flox, file = "filename.robj")

print("setting fin cluster identities")
fin.subclust$cluster.subclust.flox <- cluster.subclust.flox
Idents(fin.subclust) <- cluster.subclust.flox

fin.subclust$stage <- factor(fin.subclust$stage, levels = c("LB", "CC"))
fin.subclust$label.embryo <- factor(fin.subclust$label.embryo, levels = c("LB6",
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
subclust.umap <- data.frame(fin.subclust@reductions$ref.umap@cell.embeddings, fin.subclust$stage, fin.subclust$label.embryo, fin.subclust$cluster.whole.fin, fin.subclust$cluster.clust.fin, fin.subclust$cluster.subclust.fin, fin.subclust$cluster.whole.flox, fin.subclust$cluster.clust.flox, fin.subclust$cluster.subclust.flox)

colnames(subclust.umap) <- c("umap.1", "umap.2", "stage", "embryo", "cluster.whole.fin", "cluster.clust.fin", "cluster.subclust.fin", "cluster.whole.flox", "cluster.clust.flox", "cluster.subclust.flox")

print("saving umap coordinates")
filename.robj <- paste0(dir.name, "/projected-umap-by-stage-embryo-cluster-fin-", i, "-projected-flox-", i, "-res-", res[i], ".Robj")
save(subclust.umap, file = filename.robj)

stage.umap <- ggplot() +
geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
    labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
filename.png <- paste0(dir.name, "/projected-umap-by-stage-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = stage.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25)
stage.umap <- ggplot() +
geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~stage, ncol = 2) +
    labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())
filename.png <- paste0(dir.name, "/projected-umap-facet-stage-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = stage.umap, filename = filename.png, device = "png", width = 6.5, height = 2.25)
embryo.umap <- ggplot() +
geom_point(data = sample(subclust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
    labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
filename.png <- paste0(dir.name, "/projected-umap-by-embryo-late-bud-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25)
embryo.umap <- ggplot() +
geom_point(data = sample(subclust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~embryo, ncol = 4) +
    scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
    labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())
filename.png <- paste0(dir.name, "/projected-umap-facet-embryo-late-bud-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 6.5, height = 3.6)
embryo.umap <- ggplot() +
geom_point(data = sample(subclust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
    labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
filename.png <- paste0(dir.name, "/projected-umap-by-embryo-cardiac-crescent-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25)
embryo.umap <- ggplot() +
geom_point(data = sample(subclust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
    facet_wrap(~embryo, ncol = 4) +
    scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
    labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())
filename.png <- paste0(dir.name, "/projected-umap-facet-embryo-cardiac-crescent-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 6.5, height = 3.6)

plot.title <- paste0(i, "Nipbl FIN/+")

whole.label <- subclust.umap %>% group_by(cluster.subclust.fin) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))
cluster.umap <- ggplot() +
  geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.subclust.fin), shape = 20, stroke = 0, size = 0.5) +
  geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster.subclust.fin), size = 2.3) +
  labs(title = plot.title, x = "UMAP 1", y = "UMAP 2", color = "Nipbl FIN/+\nCluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))
cluster.umap
filename.png <- paste0(dir.name, "/projected-umap-by-fin-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")
whole.label <- subclust.umap %>% group_by(cluster.subclust.flox) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))
cluster.umap <- ggplot() +
  geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.subclust.flox), shape = 20, stroke = 0, size = 0.5) +
  geom_text(data = whole.label, mapping = aes(x = umap.1, y = umap.2, label = cluster.subclust.flox), size = 2.3) +
  labs(title = plot.title, x = "UMAP 1", y = "UMAP 2", color = "Nipbl Flox/+\nCluster") +
  guides(color = guide_legend(override.aes = list(size = 2.3))) +
  theme_classic(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))
cluster.umap
filename.png <- paste0(dir.name, "/projected-umap-by-flox-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

rm(flox.subclust)
rm(fin.anchors)
rm(fin.projected)

print("saving fin projected")
filename.robj <- paste0(dir.name, "/fin-", i, "-projected-flox-", i, "-res-", res[i], ".Robj")
save(fin.subclust, file = filename.robj)

print("preparing to find DEGs")
fin.subclust <- PrepSCTFindMarkers(fin.subclust, assay = sct.assay.name, verbose = TRUE)

print("finding DEGs")
diff.genes <- FindAllMarkers(fin.subclust, assay = sct.assay.name, logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.subclust)[2]/length(levels(fin.subclust$cluster.subclust.flox)), return.thresh = 0.05)

print("saving DEGs")
filename.robj <- paste0(dir.name, "/projected-diff-genes-fin-", i, "-projected-flox-", i, "-res-", res[i], ".Robj")
save(diff.genes, file = filename.robj)

filename.csv <- paste0(dir.name, "/projected-diff-genes-fin-", i, "-projected-flox-", i, "-res-", res[i], ".csv")
write.csv(x = diff.genes, file = filename.csv)

top.genes <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(fin.subclust$cluster.subclust.flox)))

heatmap <- DoHeatmap(fin.subclust, assay = sct.assay.name, slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
  scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
  scale_x_discrete(position = "top") +
  labs(title = plot.title, x = "Nipbl Flox/+ Cluster", y = "Top Genes", fill = "Standardized\nLog1p UMI") +
  guides(color = FALSE) +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))
filename.png <- paste0(dir.name, "/heatmap-diff-genes-fin-", i, "-projected-flox-", i, "-res-", res[i], ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

}
