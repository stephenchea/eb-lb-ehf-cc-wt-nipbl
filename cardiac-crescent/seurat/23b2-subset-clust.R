library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

for (i in c("2.2", "2.1")) {
        
    print(i)
    
    print("loading neighbors")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120")
    load("flox-2-neighbors.Robj")

    print("loading cluster identitites")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-2-res-0.056")
    load("cluster-per-cell-flox-2-res-0.056.Robj")

    print("setting cluster identities")
    flox.clust$cluster.cluster <- cluster.cluster
    Idents(flox.clust) <- cluster.cluster
    
    print("removing assays")
    DefaultAssay(flox.clust) <- "RNA"
    flox.clust[["SCT_2"]] <- NULL
    flox.clust[["INT_2"]] <- NULL
    
    print("subsetting cluster")
    flox.subclust <- subset(flox.clust, idents = i)
    
    remove(flox.clust)
    
    print("number of cells per embryo")
    print(table(flox.subclust$label.embryo))

    min.cells <- min(table(flox.subclust$label.embryo))
    print(min.cells)
    
    n <- str_replace_all(i,"[.]","_")
    sct.assay.name <- paste0("SCT_", n)
    int.assay.name <- paste0("INT_", n)

    print("applying sc transform to entire cluster")
    flox.subclust <- SCTransform(object = flox.subclust, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(flox.subclust$label.embryo))), n_genes = round(mean(flox.subclust$nFeature_RNA)), new.assay.name = sct.assay.name, verbose = TRUE)

    print("splitting by embryo")
    embryo.flox <- SplitObject(flox.subclust, split.by = "label.embryo")

    #open bracket 2
    print("applying sc transform across embryos")
    for (t in 1:length(embryo.flox)) {
        
        embryo.flox[[t]] <- SCTransform(embryo.flox[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.flox[[t]]$nFeature_RNA)), verbose = TRUE)
        
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

    for (a in 1:length(embryo.flox)) {
        embryo.flox[[a]] <- RunPCA(embryo.flox[[a]], assay = "SCT_EMBRYO", features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)
    }

    anchors.flox <- FindIntegrationAnchors(object.list = embryo.flox, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca", k.score = min(c(30, min.cells*0.95)), k.filter = min(c(200, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    flox.subclust <- IntegrateData(anchorset = anchors.flox, normalization.method = "SCT", features.to.integrate = feat.integr, new.assay.name = int.assay.name, k.weight = min(c(100, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    DefaultAssay(flox.subclust) <- int.assay.name

    print("performing principal component reduction")
    flox.subclust <- RunPCA(flox.subclust, assay = int.assay.name, features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)

    print("grabbing deviations of principal components")
    stdev.pca <- flox.subclust@reductions$pca@stdev
    dens.stdev <- density(stdev.pca)
    peak.stdev <- findpeaks(x = dens.stdev$y*-1, npeaks = 1)[1, 2]
    elbow.stdev <- dens.stdev$x[peak.stdev]
    pc.number <- sum(stdev.pca > elbow.stdev)
    print(pc.number)

    filename.robj <- paste0("sig-pc-flox-", i, ".Robj")
    save(object = pc.number, file = filename.robj)

    ###############################################
    #umap and neighbors

    flox.subclust <- RunUMAP(flox.subclust, reduction = "pca", dims = 1:pc.number, return.model = TRUE)
    filename.robj <- paste0("flox-", n, "-umap.Robj")

    flox.subclust$stage <- factor(flox.subclust$stage, levels = c("LB", "CC"))
    flox.subclust$label.embryo <- factor(flox.subclust$label.embryo, levels = c("LB1",
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

    subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust$stage, flox.subclust$label.embryo)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "stage", "embryo")
    filename.robj <- paste0("umap-by-stage-embryo-flox-", i, ".Robj")
    save(subclust.umap, file = filename.robj)

    stage.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
        
    filename.robj <- paste0("umap-by-stage-flox-", i, ".png")
    #ggsave(plot = stage.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    stage.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~stage, ncol = 2) +
        labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-stage-flox-", i, ".png")
    #ggsave(plot = stage.umap, filename = filename.robj, device = "png", width = 6.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(1, 15, 3, 17, 5)]) +
        labs(title = "Late Bud, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.robj <- paste0("umap-by-embryo-late-bud-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(1, 15, 3, 17, 5)]) +
        labs(title = "Late Bud, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-embryo-late-bud-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 6.5, height = 3.6)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(25, 13, 27, 2, 16, 4, 18, 6)]) +
        labs(title = "Cardiac Crescent, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.robj <- paste0("umap-by-embryo-cardiac-crescent-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(25, 13, 27, 2, 16, 4, 18, 6)]) +
        labs(title = "Cardiac Crescent, Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-embryo-cardiac-crescent-flox-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 6.5, height = 3.6)

    name.graph <- paste0("GRAPH_", n)
    flox.subclust <- FindNeighbors(flox.subclust, reduction = "pca", dims = 1:pc.number, graph.name = c(name.graph, name.graph))

    filename.robj <- paste0("flox-", i, "-neighbors.Robj")
    save(flox.subclust, file = filename.robj)

    ####################################
    #cluster iterations

    print("creating new directory for cluster iterations")
    direct.name <- paste0("flox-", i, "-cluster-iters")
    dir.create(direct.name)

    print("setting new working directory")
    direct.path <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-2-res-0.056/", direct.name)
    setwd(direct.path)

    name.graph <- paste0("GRAPH_", n)

    #open bracket 3
    print("begin clustering iterations")
    for (t in seq(from = 0, to = 0.5, by = 0.001)) {

    flox.subclust <- FindClusters(flox.subclust, resolution = t, graph.name = name.graph, verbose = TRUE)

    print("saving cluster identities")
    renumbered.cluster.ids <- 1:length(levels(flox.subclust))
    names(renumbered.cluster.ids) <- levels(flox.subclust)
    flox.subclust <- RenameIdents(flox.subclust, renumbered.cluster.ids)
    cluster.subclust <- flox.subclust@active.ident
    filename.robj <- paste0("cluster-per-cell-flox-", i, "-res-", sprintf(fmt = "%.3f", t), ".Robj")
    save(cluster.subclust, file = filename.robj)

    print("generating umap of clusters")
    subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust@active.ident)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "cluster")

    subclust.label <- subclust.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

    cluster.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), shape = 20, stroke = 0, size = 0.5) +
    geom_text(data = subclust.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 2.3) +
    labs(title = paste0(i, ", Nipbl Flox/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-cluster-flox-", i, "-res-", sprintf(fmt = "%.3f", t), ".png")
    #ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

    }

    remove(flox.subclust)

}
