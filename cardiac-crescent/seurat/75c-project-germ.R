library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

#begin loop across germ layers
for (i in c("MES")) {
    
    print("loading flox cells")
    filename.robj <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-", i, "-umap.Robj")
    load(file = filename.robj)

    print("loading flox pc's")
    filename.robj <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/sig-pc-flox-", i, ".Robj")
    load(file = filename.robj)

    print("setting flox cluster identities")
    flox.germ$cluster.germ <- flox.germ$cluster.whole
    flox.germ$cluster.germ <- factor(flox.germ$cluster.germ)
    Idents(flox.germ) <- flox.germ$cluster.germ
    
    print("setting working directory")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120")
    
    print("loading fin cells")
    load("fin-integr-projected-flox-integr-res-0.0120.Robj")

    print("setting fin cluster identities")
    load("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151/umap-by-cluster-fin-integr-res-0.0151-28.Robj")
    cluster.subclust.fin <- umap.whole$cluster.whole
    fin.integr$cluster.subclust.fin <- cluster.subclust.fin
    
    print("setting fin projected cluster identities")
    load("projected-umap-by-stage-embryo-cluster-fin-integr-projected-flox-integr-res-0.0120-29.Robj")
    fin.integr$cluster.whole.flox <- whole.umap$cluster.whole.flox

    print("setting fin projected germ layer identities")
    load("projected-umap-by-germ-layer-fin-integr-projected-flox-integr-res-0.0120-29.Robj")
    fin.integr$germ.whole.flox <- whole.umap$germ.whole.flox
    
    print("setting germ layer identity")
    Idents(fin.integr) <- fin.integr$germ.whole.flox
    
    print(paste0("subsetting cluster ", i))
    fin.germ <- subset(fin.integr, idents = i)
    
    remove(fin.integr)

    print("number of cells per embryo")
    print(table(fin.germ$label.embryo))

    print("smallest embryo")
    min.cells <- min(table(fin.germ$label.embryo))
    print(min.cells)

    print("generating assay names")
    sct.assay.name <- paste0("SCT_", i)
    int.assay.name <- paste0("INT_", i)

    print("performing sctransform on entire cluster")
    fin.germ <- SCTransform(object = fin.germ, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.germ$label.embryo))), n_genes = round(mean(fin.germ$nFeature_RNA)), new.assay.name = sct.assay.name, verbose = TRUE)

    print("splitting by embryo")
    embryo.fin <- SplitObject(fin.germ, split.by = "label.embryo")

    #open bracket 2
    print("performing sc transform across embryos")
    for (t in 1:length(embryo.fin)) {
        
        embryo.fin[[t]] <- SCTransform(embryo.fin[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = TRUE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.fin[[t]]$nFeature_RNA)), verbose = TRUE)
        
        DefaultAssay(embryo.fin[[t]]) <- "SCT_EMBRYO"
        
    }
    #close open bracket 2

    print("grabbing variable genes from flox")
    all.var.genes <- flox.germ@assays[[int.assay.name]]@var.features
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

    print("performing pca across all embryos")
    for (r in 1:length(embryo.fin)) {
        embryo.fin[[r]] <- RunPCA(embryo.fin[[r]], assay = "SCT_EMBRYO", features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)
    }

    print("finding integration anchors across all embryos")
    anchors.fin <- FindIntegrationAnchors(object.list = embryo.fin, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca", k.score = min(c(30, min.cells*0.95)), k.filter = min(c(200, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    print("integrating embryos")
    fin.germ <- IntegrateData(anchorset = anchors.fin, normalization.method = "SCT", features.to.integrate = feat.integr, new.assay.name = int.assay.name, k.weight = min(c(100, min.cells*0.82)), dims = 1:min(c(30, min.cells*0.95)))

    print("setting default assay to integrated assay")
    DefaultAssay(fin.germ) <- int.assay.name

    fin.germ[["SCT_EMBRYO"]] <- NULL
    
    fin.germ$cluster.germ.fin <- fin.germ$cluster.subclust.fin
    fin.germ$cluster.germ.fin <- factor(fin.germ$cluster.germ.fin)
    
    fin.germ$cluster.germ.flox <- fin.germ$cluster.whole.flox
    fin.germ$cluster.germ.flox <- factor(fin.germ$cluster.germ.flox)
    
    print("creating new directory")
    dir.name <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)))
    dir.create(dir.name)
    
    print("setting working directory")
    setwd(dir.name)

    print("finding transfer anchors")
    fin.anchors <- FindTransferAnchors(reference = flox.germ, query = fin.germ, normalization.method = "SCT", reference.assay = int.assay.name, query.assay = int.assay.name,  reduction = "pcaproject", project.query = FALSE, reference.reduction = "pca", npcs = NULL, dims = 1:pc.number, max.features = length(rownames(flox.germ)))

    print("integrating embeddings")
    fin.germ <- IntegrateEmbeddings(anchorset = fin.anchors, new.reduction.name = "ref.pca", reference = flox.germ, query = fin.germ)

    print("projecting umap")
    fin.germ <- ProjectUMAP(query = fin.germ, query.dims = 1:pc.number, reference = flox.germ, reference.dims = 1:pc.number, query.reduction = "ref.pca", reference.reduction = "pca", reduction.model = "umap")
    
    fin.germ$stage <- factor(fin.germ$stage, levels = c("LB", "CC"))
    fin.germ$label.embryo <- factor(fin.germ$label.embryo, levels = c("LB6",
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
    
    print("setting fin cluster identities")
    Idents(fin.germ) <- fin.germ$cluster.germ.flox
    
    print("saving fin germ")
    filename.robj <- filename.robj <- paste0(dir.name, "/fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".Robj")
    save(fin.germ, file = filename.robj)

    cluster.germ.flox <- fin.germ$cluster.germ.flox
    filename.robj <- paste0("cluster-per-cell-fin-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".Robj")
    save(cluster.germ.flox, file = filename.robj)
    
    print("grabbing umap coordinates")
    germ.umap <- data.frame(fin.germ@reductions$ref.umap@cell.embeddings, fin.germ$stage, fin.germ$label.embryo, fin.germ$cluster.germ.fin, fin.germ$cluster.germ.flox)

    colnames(germ.umap) <- c("umap.1", "umap.2", "stage", "embryo", "cluster.germ.fin", "cluster.germ.flox")

    print("saving umap coordinates")
    filename.robj <- paste0(dir.name, "/projected-umap-by-stage-embryo-cluster-fin-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".Robj")
    save(germ.umap, file = filename.robj)

    plot.title <- paste0(i, ", Nipbl FIN/+")
    
    stage.umap <- ggplot() +
    geom_point(data = sample(germ.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
        
            filename.png <- paste0(dir.name, "/projected-umap-by-stage-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = stage.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25)
    
    stage.umap <- ggplot() +
    geom_point(data = sample(germ.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~stage, ncol = 2) +
        labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())
            filename.png <- paste0(dir.name, "/projected-umap-facet-stage-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = stage.umap, filename = filename.png, device = "png", width = 6.5, height = 2.25)
    
    embryo.umap <- ggplot() +
    geom_point(data = sample(germ.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
        labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
            filename.png <- paste0(dir.name, "/projected-umap-by-embryo-late-bud-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25)
    
    embryo.umap <- ggplot() +
    geom_point(data = sample(germ.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
        labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())
            filename.png <- paste0(dir.name, "/projected-umap-facet-embryo-late-bud-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 6.5, height = 3.6)
    
    embryo.umap <- ggplot() +
    geom_point(data = sample(germ.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
        labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
            filename.png <- paste0(dir.name, "/projected-umap-by-embryo-cardiac-crescent-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = embryo.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25)
    
    embryo.umap <- ggplot() +
    geom_point(data = sample(germ.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
        labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())
            filename.png <- paste0(dir.name, "/projected-umap-facet-embryo-cardiac-crescent-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = embryo.umap, filename = filename.png, width = 6.5, height = 3.6)
        
    germ.label <- germ.umap %>% group_by(cluster.germ.fin) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))
    
    cluster.umap <- ggplot() +
      geom_point(data = sample(germ.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.germ.fin), shape = 20, stroke = 0, size = 0.5) +
      geom_text(data = germ.label, mapping = aes(x = umap.1, y = umap.2, label = cluster.germ.fin), size = 2.3) +
      labs(title = plot.title, x = "UMAP 1", y = "UMAP 2", color = "Nipbl FIN/+\nCluster") +
      guides(color = guide_legend(override.aes = list(size = 2.3))) +
      theme_classic(base_size = 7) +
      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))

            filename.png <- paste0(dir.name, "/projected-umap-by-fin-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")
    
    germ.label <- germ.umap %>% group_by(cluster.germ.flox) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))
    cluster.umap <- ggplot() +
      geom_point(data = sample(germ.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.germ.flox), shape = 20, stroke = 0, size = 0.5) +
      geom_text(data = germ.label, mapping = aes(x = umap.1, y = umap.2, label = cluster.germ.flox), size = 2.3) +
      labs(title = plot.title, x = "UMAP 1", y = "UMAP 2", color = "Nipbl Flox/+\nCluster") +
      guides(color = guide_legend(override.aes = list(size = 2.3))) +
      theme_classic(base_size = 7) +
      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(7, "pt"))
    cluster.umap
            filename.png <- paste0(dir.name, "/projected-umap-by-flox-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")
    
    rm(flox.germ)
    
    print("preparing to find DEGs")
    fin.germ <- PrepSCTFindMarkers(fin.germ, assay = sct.assay.name, verbose = TRUE)
    
    print("finding DEGs")
    diff.genes <- FindAllMarkers(fin.germ, assay = sct.assay.name, logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.germ)[2]/length(levels(fin.germ$cluster.germ.flox)), return.thresh = 0.05)
    
    print("saving DEGs")
    filename.robj <- paste0(dir.name, "/projected-diff-genes-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".Robj")
    save(diff.genes, file = filename.robj)
    
    filename.csv <- paste0(dir.name, "/projected-diff-genes-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".csv")
    write.csv(x = diff.genes, file = filename.csv)
    
    top.genes <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(fin.germ$cluster.germ.flox)))
        
    heatmap <- DoHeatmap(fin.germ, assay = sct.assay.name, slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
      scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
      scale_x_discrete(position = "top") +
      labs(title = plot.title, x = "Nipbl Flox/+ Cluster", y = "Top Genes", fill = "Standardized\nLog1p UMI") +
      guides(color = FALSE) +
      theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))
            filename.png <- paste0(dir.name, "/projected-heatmap-diff-genes-fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ.flox)), ".png")
    ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

}
