library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

#set working directory
setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18")

#load germ layer identities
load("projected-umap-by-germ-layer-fin-integr-projected-flox-integr-clust-5-18.Robj")

#begin loop across germ layers
for (i in c("ECT", "MES", "END")) {
    
no.clusters <- length(levels(factor(whole.umap %>% filter(germ.whole == i) %>% pull(wildtype.whole))))

print("creating new directory")
direct.name <- paste0("fin-", i, "-clust-", no.clusters)
dir.create(direct.name)

print("set new directory")
direct.path <- paste0("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18/", direct.name)
setwd(direct.path)

filename.robj <- paste0("fin-", i, "-umap.Robj")
load(file = filename.robj)

assay.name <- paste0("SCT_", i)

fin.germ <- PrepSCTFindMarkers(fin.germ, assay = assay.name, verbose = TRUE)

diff.genes <- FindAllMarkers(fin.germ, assay = assay.name, slot = "data", only.pos = TRUE, test.use = "wilcox", min.pct = 0, min.diff.pct = -Inf,  logfc.threshold = 0, return.thresh = 0.5, max.cells.per.ident = dim(fin.germ)[2]/length(levels(fin.germ$wildtype.whole)), verbose = TRUE)

filename.robj <- paste0("diff-genes-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".Robj")
save(diff.genes, file = filename.robj)

filename.csv <- paste0("diff-genes-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".csv")
write.csv(x = diff.genes, file = filename.csv)

top.diffs <- diff.genes %>% group_by(cluster) %>% slice_head(n = 79/length(levels(fin.germ$wildtype.whole)))

heatmap <- DoHeatmap(fin.germ, assay = assay.name, slot = "scale.data", features = as.vector(top.diffs$gene), group.by = "ident", group.bar = TRUE, label = TRUE, size = 2.3, angle = 0, hjust = 0.5, draw.lines = TRUE) +
scale_fill_gradientn(colors = rev(brewer.pal(n=5, name = "RdBu")), na.value = "white") +
scale_x_discrete(position = "top") +
labs(title = paste0(i, ", Nipbl FIN/+"), x = "Cluster", y = "Top Genes", fill = "Centered\nNorm UMI") +
guides(color = FALSE) +
theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7), plot.title = element_text(hjust = 0.5, size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

filename.png <- paste0("heatmap-diff-genes-fin-", i, "-clust-", length(levels(fin.germ$wildtype.whole)), ".png")
ggsave(plot = heatmap, filename = filename.png, device = "png", units = "in", width = 6.5, height = 9)

setwd("/dfs6/pub/schea2/20210505-seurat/fin-integr-projected-flox-integr-clust-5-18")

}
