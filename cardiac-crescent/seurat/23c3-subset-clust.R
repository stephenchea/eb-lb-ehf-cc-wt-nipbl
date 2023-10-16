library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

for (i in c("1.3", "1.1")) {
        
    print(i)

    print("set working directory")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-1-res-0.074")
    
    n <- str_replace_all(i,"[.]","_")

    print("loading neighbors")
    filename.robj <- paste0("flox-", i, "-neighbors.Robj")
    load(file = filename.robj)

    print("creating new directory for cluster iterations")
    direct.name <- paste0("flox-", i, "-cluster-iters")
    dir.create(direct.name)

    print("setting new working directory")
    direct.path <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/flox-1-res-0.074/", direct.name)
    setwd(direct.path)

    name.graph <- paste0("GRAPH_", n)

    print("begin clustering iterations")
    for (t in seq(from = 0.5, to = 1, by = 0.001)) {

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
