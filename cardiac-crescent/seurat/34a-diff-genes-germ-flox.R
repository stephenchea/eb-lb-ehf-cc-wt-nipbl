library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

#begin loop across germ layers
for (i in c("ECT")) {

    print(i)
    
    print("setting working directory")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120")

    print("loading umap")
    filename.robj <- paste0("flox-", i, "-umap.Robj")
    load(file = filename.robj)
    
    flox.germ$cluster.germ <- flox.germ$cluster.whole
    flox.germ$cluster.germ <- factor(flox.germ$cluster.germ)
    
    Idents(flox.germ) <- flox.germ$cluster.germ
    
    print("creating new directory")
    name.direct <- paste0("flox-", i, "-", length(levels(flox.germ$cluster.germ)))
    dir.create(name.direct)
    
    print("setting working directory")
    path.direct <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/", name.direct)
    setwd(path.direct)
    
    cluster.germ <- flox.germ$cluster.germ

    filename.robj <- paste0("cluster-per-cell-flox-", i, "-", length(levels(flox.germ$cluster.germ)), ".Robj")
    save(cluster.germ, file = filename.robj)
    
    print("saving umap coordinates across cluster")
    clust.umap <- data.frame(flox.germ@reductions$umap@cell.embeddings, flox.germ$stage, flox.germ$label.embryo, flox.germ$cluster.germ)
    colnames(clust.umap) <- c("umap.1", "umap.2", "stage", "embryo", "cluster")
    filename.robj <- paste0("umap-by-stage-embryo-cluster-flox-", i, "-", length(levels(flox.germ$cluster.germ)), ".Robj")
    save(clust.umap, file = filename.robj)
    
    cluster.umap <- ggplot() +
    geom_point(data = sample(clust.umap), mapping = aes(x = umap.1, y = umap.2, color = cluster.whole), shape = 20, stroke = 0, size = 0.5) +
        labs(title = "Nipbl Flox/+", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
        guides(color = guide_legend(override.aes = list(size = 2.3))) +
        theme_classic(base_size = 7) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.key.size = unit(7, "pt"))
        
    filename.robj <- paste0("umap-by-cluster-flox-", i, "-", length(levels(flox.germ$cluster.germ)), ".png")
    #ggsave(plot = cluster.umap, filename = filename.robj, device = "png", width = 3.5, height = 2.25)
    
    print("generating assay names")
    name.assay <- paste0("SCT_", i)

    print("preparing to find diff genes")
    flox.germ <- PrepSCTFindMarkers(flox.germ, assay = name.assay, verbose = TRUE)

    print("finding diff genes")
    diff.genes <- FindAllMarkers(flox.germ, assay = name.assay, logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(flox.germ)[2]/length(levels(flox.germ$cluster.germ)), return.thresh = 0.05)

    filename.robj <- paste0("diff-genes-flox-", i, "-", length(levels(flox.germ$cluster.germ)), ".Robj")
    save(diff.genes, file = filename.robj)

    filename.csv <- paste0("diff-genes-flox-", i, "-", length(levels(flox.germ$cluster.germ)), ".csv")
    write.csv(x = diff.genes, file = filename.csv)

    top.genes <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(flox.germ$cluster.germ)))

    heatmap <- DoHeatmap(flox.germ, assay = name.assay, slot = "scale.data", features = as.vector(top.genes$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
    scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
    scale_x_discrete(position = "top") +
    labs(title = paste0(i, ", Nipbl Flox/+"), fill = "Standardized\nLog1p UMI") +
    guides(color = FALSE) +
    theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

    filename.png <- paste0("heatmap-diff-genes-flox-", i, "-", length(levels(flox.germ$cluster.germ)), ".png")
    #ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

}
