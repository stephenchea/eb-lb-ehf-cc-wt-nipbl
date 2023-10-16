library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

###########W1
setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5")
load("flox-W1-neighbors.Robj")

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W1-clust-2")
load("cluster-per-cell-flox-W1-clust-2.Robj")

flox.clust$wildtype.clust <- flox.clust.active.ident
Idents(flox.clust) <- flox.clust.active.ident

save(flox.clust, file = "flox-W1-clust-2.Robj")

for (i in c("W1.1")) {
        
    print(i)
    print("subsetting cluster")
    flox.subclust <- subset(flox.clust, idents = i)
    
    n <- str_replace_all(i,"[.]","_")
    assay.name <- paste0("SCT_", n)

    print("applying sc transform to entire cluster")
    flox.subclust <- SCTransform(object = flox.subclust, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = round(dim(flox.subclust)[2]*.5), n_genes = round(dim(flox.subclust)[1]*.1), new.assay.name = assay.name, verbose = TRUE)

    print("splitting by embryo")
    flox.embryo <- SplitObject(flox.subclust, split.by = "embryo.id")

    #open bracket 2
    print("applying sc transform across embryos")
    for (t in 1:length(flox.embryo)) {
        
        flox.embryo[[t]] <- SCTransform(flox.embryo[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = NULL, n_genes = round(dim(flox.embryo[[t]])[1]*.1), verbose = TRUE)
        
        DefaultAssay(flox.embryo[[t]]) <- "SCT_EMBRYO"
        
    }
    #close open bracket 2

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

    print("performing principal component reduction")
    flox.subclust <- RunPCA(flox.subclust, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

    print("grabbing deviations of principal components")
    stdev.pca <- flox.clust@reductions$pca@stdev
    var.exp <- stdev.pca^2/sum(stdev.pca^2)

    print("set number of principal components")
    d.var.exp <- var.exp - median(var.exp)
    pc.number <- sum(d.var.exp > mad(var.exp)*2)

    filename.robj <- paste0("sig-pc-flox-", n, ".Robj")
    save(object = pc.number, file = filename.robj)

    ###############################################
    #umap and neighbors

    flox.subclust <- RunUMAP(flox.subclust, reduction = "pca", dims = 1:pc.number, return.model = TRUE)
    filename.robj <- paste0("flox-", n, "-umap.Robj")

    subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust$embryo.id)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "embryo.id")
    filename.robj <- paste0("umap-by-embryo-flox-", n, "-umap.Robj")
    save(subclust.umap, file = filename.robj)

    umap.embryo <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = embryo.id), size = 0.1) +
    scale_color_manual(labels = c("EP1", "EP2", "EP3", "EP4", "EP5"), values = hue_pal()(12)[c(1, 7, 3, 9, 5)]) +
    labs(title = paste0(i, ", Nipbl Flox/+"), x = "UMAP 1", y = "UMAP 2", color = "Embryo\nPair") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-embryo-flox-", n, ".png")
    ggsave(plot = umap.embryo, filename = filename.png, device = "png", width = 6.5, height = 4.5)

    name.graph <- paste0("GRAPH_", n)
    flox.subclust <- FindNeighbors(flox.subclust, reduction = "pca", dims = 1:pc.number, graph.name = c(name.graph, name.graph))

    filename.robj <- paste0("flox-", n, "-neighbors.Robj")
    save(flox.subclust, file = filename.robj)


    ####################################
    #cluster iterations

    print("creating new directory for cluster iterations")
    direct.name <- paste0("flox-", n, "-cluster-iters")
    dir.create(direct.name)

    print("setting new working directory")
    direct.path <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W1-clust-2/", direct.name)
    setwd(direct.path)

    name.graph <- paste0("GRAPH_", n)

    res = 0.1

    #open bracket 3
    print("begin clustering iterations")
    for (t in 2:10) {
        
    flox.subclust <- FindClusters(flox.subclust, resolution = res, graph.name = name.graph, random.seed = 1, verbose = TRUE)

    cluster.number = t
    step = 0.1
    loop.count = 1
    loop.vector = vector()

    while (max(as.numeric(flox.subclust@active.ident[])) != cluster.number & (loop.count < 1000))
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
        if(max(as.numeric(flox.subclust@active.ident[])) < cluster.number)
        {
            res = res + step
            flox.subclust <- FindClusters(object = flox.subclust, resolution = res, graph.name = name.graph, random.seed = 1, verbose = TRUE)
            loop.vector <- c(loop.vector, max(as.numeric(flox.subclust@active.ident[])))
        }
        
        #If the max cluster is GREATER than 15, we need a slightly lower resolution
        if(max(as.numeric(flox.subclust@active.ident[])) > cluster.number)
        {
            res = res - step
            flox.subclust <- FindClusters(object = flox.subclust, resolution = res, graph.name = name.graph,  random.seed = 1, verbose = TRUE)
            loop.vector <- c(loop.vector, max(as.numeric(flox.subclust@active.ident[])))
        }
        
        print(paste0("res = ", res))
    }

    print(res)
    print(length(levels(flox.subclust)))

    print("saving cluster identities")
    renumbered.cluster.ids <- 1:t
    names(renumbered.cluster.ids) <- levels(flox.subclust)
    flox.subclust <- RenameIdents(flox.subclust, renumbered.cluster.ids)
    flox.subclust.active.ident <- flox.subclust@active.ident
    filename.robj <- paste0("cluster-per-cell-flox-", n, "-clust-", t, ".Robj")
    save(flox.subclust.active.ident, file = filename.robj)

    print("generating umap of clusters")
    subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust@active.ident)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "cluster")

    subclust.label <- subclust.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

    cluster.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), size = 0.1) +
    geom_text(data = subclust.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 3.5) +
    labs(title = paste0(i, ", Nipbl Flox/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-cluster-flox-", n, "-clust-", t, ".png")
    ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

    res = res + step
    }

    setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W1-clust-2")

}

###########W2
setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5")
load("flox-W2-neighbors.Robj")

setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-3")
load("cluster-per-cell-flox-W2-clust-3.Robj")

flox.clust$wildtype.clust <- flox.clust.active.ident
Idents(flox.clust) <- flox.clust.active.ident

save(flox.clust, file = "flox-W2-clust-3.Robj")

for (i in c("W2.1", "W2.2", "W2.3")) {
        
    print(i)
    print("subsetting cluster")
    flox.subclust <- subset(flox.clust, idents = i)
    
    n <- str_replace_all(i,"[.]","_")
    assay.name <- paste0("SCT_", n)

    print("applying sc transform to entire cluster")
    flox.subclust <- SCTransform(object = flox.subclust, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = round(dim(flox.subclust)[2]*.5), n_genes = round(dim(flox.subclust)[1]*.1), new.assay.name = assay.name, verbose = TRUE)

    print("splitting by embryo")
    flox.embryo <- SplitObject(flox.subclust, split.by = "embryo.id")

    #open bracket 2
    print("applying sc transform across embryos")
    for (t in 1:length(flox.embryo)) {
        
        flox.embryo[[t]] <- SCTransform(flox.embryo[[t]], assay = "RNA", new.assay.name = "SCT_EMBRYO", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = NULL, n_genes = round(dim(flox.embryo[[t]])[1]*.1), verbose = TRUE)
        
        DefaultAssay(flox.embryo[[t]]) <- "SCT_EMBRYO"
        
    }
    #close open bracket 2

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

    print("performing principal component reduction")
    flox.subclust <- RunPCA(flox.subclust, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

    print("grabbing deviations of principal components")
    stdev.pca <- flox.clust@reductions$pca@stdev
    var.exp <- stdev.pca^2/sum(stdev.pca^2)

    print("set number of principal components")
    d.var.exp <- var.exp - median(var.exp)
    pc.number <- sum(d.var.exp > mad(var.exp)*2)

    filename.robj <- paste0("sig-pc-flox-", n, ".Robj")
    save(object = pc.number, file = filename.robj)

    ###############################################
    #umap and neighbors

    flox.subclust <- RunUMAP(flox.subclust, reduction = "pca", dims = 1:pc.number, return.model = TRUE)
    filename.robj <- paste0("flox-", n, "-umap.Robj")

    subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust$embryo.id)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "embryo.id")
    filename.robj <- paste0("umap-by-embryo-flox-", n, "-umap.Robj")
    save(subclust.umap, file = filename.robj)

    umap.embryo <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = embryo.id), size = 0.1) +
    scale_color_manual(labels = c("EP1", "EP2", "EP3", "EP4", "EP5"), values = hue_pal()(12)[c(1, 7, 3, 9, 5)]) +
    labs(title = paste0(i, ", Nipbl Flox/+"), x = "UMAP 1", y = "UMAP 2", color = "Embryo\nPair") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-embryo-flox-", n, ".png")
    ggsave(plot = umap.embryo, filename = filename.png, device = "png", width = 6.5, height = 4.5)

    name.graph <- paste0("GRAPH_", n)
    flox.subclust <- FindNeighbors(flox.subclust, reduction = "pca", dims = 1:pc.number, graph.name = c(name.graph, name.graph))

    filename.robj <- paste0("flox-", n, "-neighbors.Robj")
    save(flox.subclust, file = filename.robj)


    ####################################
    #cluster iterations

    print("creating new directory for cluster iterations")
    direct.name <- paste0("flox-", n, "-cluster-iters")
    dir.create(direct.name)

    print("setting new working directory")
    direct.path <- paste0("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-3/", direct.name)
    setwd(direct.path)

    name.graph <- paste0("GRAPH_", n)

    res = 0.1

    #open bracket 3
    print("begin clustering iterations")
    for (t in 2:10) {
        
    flox.subclust <- FindClusters(flox.subclust, resolution = res, graph.name = name.graph, random.seed = 1, verbose = TRUE)

    cluster.number = t
    step = 0.1
    loop.count = 1
    loop.vector = vector()

    while (max(as.numeric(flox.subclust@active.ident[])) != cluster.number & (loop.count < 1000))
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
        if(max(as.numeric(flox.subclust@active.ident[])) < cluster.number)
        {
            res = res + step
            flox.subclust <- FindClusters(object = flox.subclust, resolution = res, graph.name = name.graph, random.seed = 1, verbose = TRUE)
            loop.vector <- c(loop.vector, max(as.numeric(flox.subclust@active.ident[])))
        }
        
        #If the max cluster is GREATER than 15, we need a slightly lower resolution
        if(max(as.numeric(flox.subclust@active.ident[])) > cluster.number)
        {
            res = res - step
            flox.subclust <- FindClusters(object = flox.subclust, resolution = res, graph.name = name.graph,  random.seed = 1, verbose = TRUE)
            loop.vector <- c(loop.vector, max(as.numeric(flox.subclust@active.ident[])))
        }
        
        print(paste0("res = ", res))
    }

    print(res)
    print(length(levels(flox.subclust)))

    print("saving cluster identities")
    renumbered.cluster.ids <- 1:t
    names(renumbered.cluster.ids) <- levels(flox.subclust)
    flox.subclust <- RenameIdents(flox.subclust, renumbered.cluster.ids)
    flox.subclust.active.ident <- flox.subclust@active.ident
    filename.robj <- paste0("cluster-per-cell-flox-", n, "-clust-", t, ".Robj")
    save(flox.subclust.active.ident, file = filename.robj)

    print("generating umap of clusters")
    subclust.umap <- data.frame(flox.subclust@reductions$umap@cell.embeddings, flox.subclust@active.ident)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "cluster")

    subclust.label <- subclust.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

    cluster.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), size = 0.1) +
    geom_text(data = subclust.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 3.5) +
    labs(title = paste0(i, ", Nipbl Flox/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-cluster-flox-", n, "-clust-", t, ".png")
    ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

    res = res + step
    }

    setwd("/dfs6/pub/schea2/20210505-seurat/flox-integr-clust-5/flox-W2-clust-3")

}
