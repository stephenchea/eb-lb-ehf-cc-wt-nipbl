library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

for (i in c("2.2", "2.1")) {
        
    print(i)
    
    print("loading neighbors")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151")
    load("fin-2-neighbors.Robj")

    print("loading cluster identitites")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151/fin-2-res-0.042")
    load("cluster-per-cell-fin-2-res-0.042.Robj")

    print("setting cluster identities")
    fin.clust$cluster.cluster <- cluster.cluster
    Idents(fin.clust) <- cluster.cluster
    
    print("removing assays")
    DefaultAssay(fin.clust) <- "RNA"
    fin.clust[["SCT_2"]] <- NULL
    fin.clust[["INT_2"]] <- NULL
    
    print("subsetting cluster")
    fin.subclust <- subset(fin.clust, idents = i)
    
    remove(fin.clust)
    
    print("number of cells per embryo")
    print(table(fin.subclust$label.embryo))

    min.cells <- min(table(fin.subclust$label.embryo))
    print(min.cells)
    
    n <- str_replace_all(i,"[.]","_")
    sct.assay.name <- paste0("SCT_", n)
    int.assay.name <- paste0("INT_", n)

    print("applying sc transform to entire cluster")
    fin.subclust <- SCTransform(object = fin.subclust, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.subclust$label.embryo))), n_genes = round(mean(fin.subclust$nFeature_RNA)), new.assay.name = sct.assay.name, verbose = TRUE)

    print("splitting by embryo")
    embryo.fin <- SplitObject(fin.subclust, split.by = "label.embryo")

    #open bracket 2
    print("applying sc transform across embryos")
    for (t in 1:length(embryo.fin)) {
        
        embryo.fin[[t]] <- SCTransform(embryo.fin[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = NULL, n_genes = round(mean(embryo.fin[[t]]$nFeature_RNA)), verbose = TRUE)
        
        DefaultAssay(embryo.fin[[t]]) <- "SCT_EMBRYO"
        
    }
    #close open bracket 2

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

    for (a in 1:length(embryo.fin)) {
        embryo.fin[[a]] <- RunPCA(embryo.fin[[a]], assay = "SCT_EMBRYO", features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)
    }

    anchors.fin <- FindIntegrationAnchors(object.list = embryo.fin, normalization.method = "SCT", anchor.features = feat.integr, reduction = "rpca", k.score = min(c(30, min.cells*0.95)), k.filter = min(c(200, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    fin.subclust <- IntegrateData(anchorset = anchors.fin, normalization.method = "SCT", features.to.integrate = feat.integr, new.assay.name = int.assay.name, k.weight = min(c(100, min.cells*0.95)), dims = 1:min(c(30, min.cells*0.95)))

    DefaultAssay(fin.subclust) <- int.assay.name

    print("performing principal component reduction")
    fin.subclust <- RunPCA(fin.subclust, assay = int.assay.name, features = feat.integr, npcs = min(c(50, min.cells*0.95)), verbose = TRUE)

    print("grabbing deviations of principal components")
    stdev.pca <- fin.subclust@reductions$pca@stdev
    dens.stdev <- density(stdev.pca)
    peak.stdev <- findpeaks(x = dens.stdev$y*-1, npeaks = 1)[1, 2]
    elbow.stdev <- dens.stdev$x[peak.stdev]
    pc.number <- sum(stdev.pca > elbow.stdev)
    
    pc.number <- max(c(pc.number, 3))
    print(pc.number)

    filename.robj <- paste0("sig-pc-fin-", i, ".Robj")
    save(object = pc.number, file = filename.robj)

    ###############################################
    #umap and neighbors

    fin.subclust <- RunUMAP(fin.subclust, reduction = "pca", dims = 1:pc.number, return.model = TRUE)
    filename.robj <- paste0("fin-", n, "-umap.Robj")

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

    subclust.umap <- data.frame(fin.subclust@reductions$umap@cell.embeddings, fin.subclust$stage, fin.subclust$label.embryo)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "stage", "embryo")
    filename.robj <- paste0("umap-by-stage-embryo-fin-", i, ".Robj")
    save(subclust.umap, file = filename.robj)

    stage.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
        
    filename.robj <- paste0("umap-by-stage-fin-", i, ".png")
    #ggsave(plot = stage.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    stage.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = stage), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~stage, ncol = 2) +
        labs(title = "Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Stage") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-stage-fin-", i, ".png")
    #ggsave(plot = stage.umap, filename = filename.robj, device = "png", width = 6.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
        labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.robj <- paste0("umap-by-embryo-late-bud-fin-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "LB")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(19, 7, 21, 9, 23, 11)]) +
        labs(title = "Late Bud, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-embryo-late-bud-fin-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 6.5, height = 3.6)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
        labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.robj <- paste0("umap-by-embryo-cardiac-crescent-fin-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)

    embryo.umap <- ggplot() +
    geom_point(data = sample(subclust.umap %>% filter(stage == "CC")), mapping = aes(x = umap.1, y = umap.2, color = embryo), shape = 20, stroke = 0, size = 0.5) +
        facet_wrap(~embryo, ncol = 4) +
        scale_color_manual(values = hue_pal()(28)[c(20, 8, 22, 10, 24, 12, 26, 14)]) +
        labs(title = "Cardiac Crescent, Nipbl FIN/+", x = "UMAP 1", y = "UMAP 2", color = "Embryo") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"), strip.background = element_blank())

    filename.robj <- paste0("umap-facet-embryo-cardiac-crescent-fin-", i, ".png")
    #ggsave(plot = embryo.umap, filename = filename.robj, device = "png", width = 6.5, height = 3.6)

    name.graph <- paste0("GRAPH_", n)
    fin.subclust <- FindNeighbors(fin.subclust, reduction = "pca", dims = 1:pc.number, graph.name = c(name.graph, name.graph))

    filename.robj <- paste0("fin-", i, "-neighbors.Robj")
    save(fin.subclust, file = filename.robj)

    ####################################
    #cluster iterations

    print("creating new directory for cluster iterations")
    direct.name <- paste0("fin-", i, "-cluster-iters")
    dir.create(direct.name)

    print("setting new working directory")
    direct.path <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-res-0.0151/fin-2-res-0.042/", direct.name)
    setwd(direct.path)

    name.graph <- paste0("GRAPH_", n)

    #open bracket 3
    print("begin clustering iterations")
    for (t in seq(from = 0, to = 1, by = 0.001)) {

    fin.subclust <- FindClusters(fin.subclust, resolution = t, graph.name = name.graph, verbose = TRUE)

    print("saving cluster identities")
    renumbered.cluster.ids <- 1:length(levels(fin.subclust))
    names(renumbered.cluster.ids) <- levels(fin.subclust)
    fin.subclust <- RenameIdents(fin.subclust, renumbered.cluster.ids)
    cluster.subclust <- fin.subclust@active.ident
    filename.robj <- paste0("cluster-per-cell-fin-", i, "-res-", sprintf(fmt = "%.3f", t), ".Robj")
    save(cluster.subclust, file = filename.robj)

    print("generating umap of clusters")
    subclust.umap <- data.frame(fin.subclust@reductions$umap@cell.embeddings, fin.subclust@active.ident)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "cluster")

    subclust.label <- subclust.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

    cluster.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), shape = 20, stroke = 0, size = 0.5) +
    geom_text(data = subclust.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 2.3) +
    labs(title = paste0(i, ", Nipbl FIN/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-cluster-fin-", i, "-res-", sprintf(fmt = "%.3f", t), ".png")
    #ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 3.5, height = 2.25, units = "in")

    }

    remove(fin.subclust)

}
