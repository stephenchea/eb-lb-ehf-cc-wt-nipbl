library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

###########M2
setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6")
load("fin-M2-neighbors.Robj")

setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6/fin-M2-clust-3")
load("cluster-per-cell-fin-M2-clust-3.Robj")

fin.clust$mutant.clust <- fin.clust.active.ident
Idents(fin.clust) <- fin.clust.active.ident

save(fin.clust, file = "fin-M2-clust-3.Robj")

for (i in c("M2.1", "M2.2")) {
        
    print(i)
    print("subsetting cluster")
    fin.subclust <- subset(fin.clust, idents = i)
    
    n <- str_replace_all(i,"[.]","_")
    assay.name <- paste0("SCT_", n)

    print("applying sc transform to entire cluster")
    fin.subclust <- SCTransform(object = fin.subclust, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 3, n_cells = round(dim(fin.subclust)[2]*.5), n_genes = round(dim(fin.subclust)[1]*.1), new.assay.name = assay.name, verbose = TRUE)

    print("splitting by embryo")
    fin.embryo <- SplitObject(fin.subclust, split.by = "embryo.id")

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
    fin.subclust <- RunPCA(fin.subclust, assay = "integrated", features = feat.integr, npcs = 50, verbose = TRUE)

    print("grabbing deviations of principal components")
    stdev.pca <- fin.clust@reductions$pca@stdev
    var.exp <- stdev.pca^2/sum(stdev.pca^2)

    print("set number of principal components")
    d.var.exp <- var.exp - median(var.exp)
    pc.number <- sum(d.var.exp > mad(var.exp)*2)

    filename.robj <- paste0("sig-pc-fin-", n, ".Robj")
    save(object = pc.number, file = filename.robj)

    ###############################################
    #umap and neighbors

    fin.subclust <- RunUMAP(fin.subclust, reduction = "pca", dims = 1:pc.number)
    filename.robj <- paste0("fin-", n, "-umap.Robj")

    subclust.umap <- data.frame(fin.subclust@reductions$umap@cell.embeddings, fin.subclust$embryo.id)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "embryo.id")
    filename.robj <- paste0("umap-by-embryo-fin-", n, "-umap.Robj")
    save(subclust.umap, file = filename.robj)

    umap.embryo <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = embryo.id), size = 0.1) +
    scale_color_manual(labels = c("EP6", "EP7", "EP8", "EP9", "EP10", "EP11"), values = hue_pal()(12)[c(2, 8, 4, 10, 6, 12)]) +
    labs(title = paste0(i, ", Nipbl FIN/+"), x = "UMAP 1", y = "UMAP 2", color = "Embryo\nPair") +
    guides(color = guide_legend(override.aes = list(size = 2.3))) +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-embryo-fin-", n, ".png")
    ggsave(plot = umap.embryo, filename = filename.png, device = "png", width = 6.5, height = 4.5)

    name.graph <- paste0("GRAPH_", n)
    fin.subclust <- FindNeighbors(fin.subclust, reduction = "pca", dims = 1:pc.number, graph.name = c(name.graph, name.graph))

    filename.robj <- paste0("fin-", n, "-neighbors.Robj")
    save(fin.subclust, file = filename.robj)


    ####################################
    #cluster iterations

    print("creating new directory for cluster iterations")
    direct.name <- paste0("fin-", n, "-cluster-iters")
    dir.create(direct.name)

    print("setting new working directory")
    direct.path <- paste0("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6/fin-M2-clust-3/", direct.name)
    setwd(direct.path)

    name.graph <- paste0("GRAPH_", n)

    res = 0.1

    #open bracket 3
    print("begin clustering iterations")
    for (t in 2:10) {
        
    fin.subclust <- FindClusters(fin.subclust, resolution = res, graph.name = name.graph, random.seed = 1, verbose = TRUE)

    cluster.number = t
    step = 0.1
    loop.count = 1
    loop.vector = vector()

    while (max(as.numeric(fin.subclust@active.ident[])) != cluster.number & (loop.count < 1000))
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
        if(max(as.numeric(fin.subclust@active.ident[])) < cluster.number)
        {
            res = res + step
            fin.subclust <- FindClusters(object = fin.subclust, resolution = res, graph.name = name.graph, random.seed = 1, verbose = TRUE)
            loop.vector <- c(loop.vector, max(as.numeric(fin.subclust@active.ident[])))
        }
        
        #If the max cluster is GREATER than 15, we need a slightly lower resolution
        if(max(as.numeric(fin.subclust@active.ident[])) > cluster.number)
        {
            res = res - step
            fin.subclust <- FindClusters(object = fin.subclust, resolution = res, graph.name = name.graph,  random.seed = 1, verbose = TRUE)
            loop.vector <- c(loop.vector, max(as.numeric(fin.subclust@active.ident[])))
        }
        
        print(paste0("res = ", res))
    }

    print(res)
    print(length(levels(fin.subclust)))

    print("saving cluster identities")
    renumbered.cluster.ids <- 1:t
    names(renumbered.cluster.ids) <- levels(fin.subclust)
    fin.subclust <- RenameIdents(fin.subclust, renumbered.cluster.ids)
    fin.subclust.active.ident <- fin.subclust@active.ident
    filename.robj <- paste0("cluster-per-cell-fin-", n, "-clust-", t, ".Robj")
    save(fin.subclust.active.ident, file = filename.robj)

    print("generating umap of clusters")
    subclust.umap <- data.frame(fin.subclust@reductions$umap@cell.embeddings, fin.subclust@active.ident)
    colnames(subclust.umap) <- c("umap.1", "umap.2", "cluster")

    umap.label <- subclust.umap %>% group_by(cluster) %>% summarize(umap.1 = median(umap.1), umap.2 = median(umap.2))

    cluster.umap <- ggplot() +
    geom_point(data = sample(subclust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster), size = 0.1) +
    geom_text(data = umap.label, mapping = aes(x = umap.1, y = umap.2, label = cluster), size = 3.5) +
    labs(title = paste0(i, ", Nipbl FIN/+"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
    theme_classic(base_size = 7) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("umap-by-cluster-fin-", n, "-clust-", t, ".png")
    ggsave(plot = cluster.umap, filename = filename.png, device = "png", width = 6.5, height = 4.5, units = "in")

    res = res + step
    }

    setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-clust-6/fin-M2-clust-3")

}

