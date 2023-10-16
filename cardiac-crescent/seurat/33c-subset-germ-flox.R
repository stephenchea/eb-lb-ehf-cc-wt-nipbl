library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

#begin loop across germ layers
for (i in c("MES")) {
    
    #set working directory
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

    #load seurat object
    load("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-neighbors.Robj")

    #set working directory
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120")

    #load cluster identities
    load("umap-by-cluster-flox-integr-res-0.0120-29.Robj")
    flox.integr$cluster.whole <- umap.whole$cluster.whole

    #load germ layer identities
    load("umap-by-germ-layer-flox-integr-res-0.0120-29.Robj")
    flox.integr$germ.whole <- umap.whole$germ.whole

    #set identity to germ layer
    Idents(flox.integr) <- umap.whole$germ.whole
    
    print(paste0("subsetting cluster ", i))
    flox.germ <- subset(flox.integr, idents = i)
    
    remove(flox.integr)

    print("number of cells per embryo")
    print(table(flox.germ$label.embryo))

    print("smallest embryo")
    min.cells <- min(table(flox.germ$label.embryo))
    print(min.cells)

    print("generating assay names")
    sct.assay.name <- paste0("SCT_", i)
    int.assay.name <- paste0("INT_", i)

    print("performing sctransform on entire cluster")
    flox.germ <- SCTransform(object = flox.germ, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(flox.germ$label.embryo))), n_genes = round(mean(flox.germ$nFeature_RNA)), new.assay.name = sct.assay.name, verbose = TRUE)

    print("splitting by embryo")
    embryo.flox <- SplitObject(flox.germ, split.by = "label.embryo")

    #open bracket 2
    print("performing sc transform across embryos")
    for (t in 1:length(embryo.flox)) {
        
        embryo.flox[[t]] <- SCTransform(embryo.flox[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = TRUE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.flox[[t]]$nFeature_RNA)), verbose = TRUE)
        
        DefaultAssay(embryo.flox[[t]]) <- "SCT_EMBRYO"
        
    }
    #close open bracket 2

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

    print("calculating number of variable genes per embryo")
    lapply(var.genes, FUN = length)

    print("finding unique variable genes among late bud embryos")
    lb.uniq.genes <- unique(c(var.genes[[1]],
    var.genes[[2]],
    var.genes[[3]],
    var.genes[[4]],
    var.genes[[5]]
    ))

    print("calculating number of unique variable genes in late bud embyros")
    length(lb.uniq.genes)

    print("counting how many times each gene occurs among late bud embryos")
    lb.count.gene <- list()
    for (e in lb.uniq.genes) {
        for(r in names(var.genes)[1:5]) {
            lb.count.gene[[e]][[r]]<- e %in% var.genes[[r]]
        }
    }

    print("ranking unique variable genes among late bud embryos")
    lb.count.rank <- lapply(X = lb.count.gene, FUN = unlist)
    lb.count.rank <- lapply(X = lb.count.rank, FUN = sum)
    lb.count.rank <- sort(unlist(lb.count.rank), decreasing = TRUE)

    print("table of rankings among late bud embryos")
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

    print("calculating number of unique variable genes in cardiac crescent embyros")
    length(cc.uniq.genes)

    print("counting how many times each gene occurs among cardiac crescent embyros")
    cc.count.gene <- list()
    for (e in cc.uniq.genes) {
        for(r in names(var.genes)[6:13]) {
            cc.count.gene[[e]][[r]]<- e %in% var.genes[[r]]
        }
    }

    print("ranking unique variable genes among cardiac crescent embyros")
    cc.count.rank <- lapply(X = cc.count.gene, FUN = unlist)
    cc.count.rank <- lapply(X = cc.count.rank, FUN = sum)
    cc.count.rank <- sort(unlist(cc.count.rank), decreasing = TRUE)

    print("table of rankings among cardiac crescent embyros")
    print(table(cc.count.rank))

    print("grabbing variable genes common to all cardiac crescent embryos")
    cc.var.genes <- names(cc.count.rank[cc.count.rank >= 8])

    print("grabbing variable genes common to late bud and cardiac crescent embryos")
    all.var.genes <- unique(c(lb.var.genes, cc.var.genes))

    print("removing mitochondrial genes from variable genes")
    all.var.genes <- str_subset(all.var.genes, "mt-", negate = TRUE)

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

    print("keeping only variable genes among common sct genes")
    feat.integr <- intersect(all.var.genes, Reduce(intersect, sct.genes))
    length(feat.integr)

    print("preparing sctransform assay for integration")
    embryo.flox <- PrepSCTIntegration(object.list = embryo.flox, anchor.features = feat.integr)

    print("performing pca across all embryos")
    for (r in 1:length(embryo.flox)) {
        embryo.flox[[r]] <- RunPCA(embryo.flox[[r]], assay = "SCT_EMBRYO", features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)
    }

    print("finding integration anchors across all embryos")
    anchors.flox <- FindIntegrationAnchors(object.list = embryo.flox, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca", k.score = min(c(30, min.cells*0.95)), k.filter = min(c(200, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    print("integrating embryos")
    flox.germ <- IntegrateData(anchorset = anchors.flox, normalization.method = "SCT", features.to.integrate = feat.integr, new.assay.name = int.assay.name, k.weight = min(c(100, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    print("setting default assay to integrated assay")
    DefaultAssay(flox.germ) <- int.assay.name

    flox.germ[["SCT_EMBRYO"]] <- NULL

    print("performing pca using integrated assay")
    flox.germ <- RunPCA(flox.germ, assay = int.assay.name, features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)

    print("grabbing deviations of principal components")
    stdev.pca <- flox.germ@reductions$pca@stdev
    dens.stdev <- density(stdev.pca)
    peak.stdev <- findpeaks(x = dens.stdev$y*-1, npeaks = 1)[1, 2]
    elbow.stdev <- dens.stdev$x[peak.stdev]
    pc.number <- sum(stdev.pca > elbow.stdev)
    pc.number <- max(c(pc.number, 3))
    print("number of principal components")
    print(pc.number)

    print("saving number of principal components")
    filename.robj <- paste0("sig-pc-flox-", i, ".Robj")
    save(object = pc.number, file = filename.robj)

    print("running umap on integrated embryos")
    flox.germ <- RunUMAP(flox.germ, reduction = "pca", dims = 1:pc.number, return.model = TRUE)

    flox.germ$stage <- factor(flox.germ$stage, levels = c("LB", "CC"))
    flox.germ$label.embryo <- factor(flox.germ$label.embryo, levels = c("LB1",
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

    print("saving umap")
    filename.robj <- paste0("flox-", i, "-umap.Robj")
    save(flox.germ, file = filename.robj)
    
    print("saving umap coordinates across stage and embryo")

    clust.umap <- data.frame(flox.germ@reductions$umap@cell.embeddings, flox.germ$stage, flox.germ$label.embryo)
    colnames(clust.umap) <- c("umap.1", "umap.2", "stage", "embryo")
    filename.robj <- paste0("umap-by-stage-embryo-flox-", i, ".Robj")
    save(clust.umap, file = filename.robj)

    print("generating umap plots by stage and embryo")

    stage.umap <- ggplot() +
    geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
        
    filename.robj <- paste0("umap-by-stage-flox-", i, ".png")
    #ggsave(plot = stage.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    stage.umap <- ggplot() +
    geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~stage, ncol = 2) +
        labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-stage-flox-", i, ".png")
    #ggsave(plot = stage.umap, filename = filename.robj, device = "png", width = 6.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(clust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(1, 15, 3, 17, 5)]) +
        labs(title = "Late Bud, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.robj <- paste0("umap-by-embryo-late-bud-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(clust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(1, 15, 3, 17, 5)]) +
        labs(title = "Late Bud, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-embryo-late-bud-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 6.5, height = 3.6)

    embryo.umap <- ggplot() +
    geom_point(data = sample(clust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(25, 13, 27, 2, 16, 4, 18, 6)]) +
        labs(title = "Cardiac Crescent, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.robj <- paste0("umap-by-embryo-cardiac-crescent-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(clust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(25, 13, 27, 2, 16, 4, 18, 6)]) +
        labs(title = "Cardiac Crescent, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-embryo-cardiac-crescent-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 6.5, height = 3.6)

}
