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
        
    print("creating new directory")
    name.direct <- paste0("flox-", i, "-", length(levels(flox.germ$cluster.germ)))
    dir.create(name.direct)
    
    print("setting working directory")
    path.direct <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/", name.direct)
    setwd(path.direct)
    
    Idents(flox.germ) <- flox.germ$stage

    #subsetting cluster stage
    flox.stage <- subset(flox.germ, idents = "LB")
    
    rm(flox.germ)
    
    print("performing sctransform on entire cluster")
    flox.stage <- SCTransform(object = flox.stage, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(flox.stage$label.embryo))), n_genes = round(mean(flox.stage$nFeature_RNA)), new.assay.name = "SCT_GERM", verbose = TRUE)
    
    Idents(flox.stage) <- flox.stage$cluster.germ

    #print("preparing to find diff genes")
    flox.stage <- PrepSCTFindMarkers(flox.stage, assay = "SCT_GERM", verbose = TRUE)

    print("finding diff genes")
    diff.genes <- FindAllMarkers(flox.stage, assay = "SCT_GERM", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(flox.stage)[2]/length(levels(flox.stage$cluster.germ)), return.thresh = 0.05)

    filename.robj <- paste0("diff-genes-LB-flox-", i, "-", length(levels(flox.stage$cluster.germ)), ".Robj")
    save(diff.genes, file = filename.robj)

    filename.csv <- paste0("diff-genes-LB-flox-", i, "-", length(levels(flox.stage$cluster.germ)), ".csv")
    write.csv(x = diff.genes, file = filename.csv)
    
    rm(flox.stage)
    
    print(i)
    
    print("setting working directory")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120")

    print("loading umap")
    filename.robj <- paste0("flox-", i, "-umap.Robj")
    load(file = filename.robj)
    
    flox.germ$cluster.germ <- flox.germ$cluster.whole
    flox.germ$cluster.germ <- factor(flox.germ$cluster.germ)
        
    print("creating new directory")
    name.direct <- paste0("flox-", i, "-", length(levels(flox.germ$cluster.germ)))
    dir.create(name.direct)
    
    print("setting working directory")
    path.direct <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120/", name.direct)
    setwd(path.direct)
    
    Idents(flox.germ) <- flox.germ$stage

    #subsetting cluster stage
    flox.stage <- subset(flox.germ, idents = "CC")
    
    rm(flox.germ)
    
    print("performing sctransform on entire cluster")
    flox.stage <- SCTransform(object = flox.stage, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(flox.stage$label.embryo))), n_genes = round(mean(flox.stage$nFeature_RNA)), new.assay.name = "SCT_GERM", verbose = TRUE)
    
    Idents(flox.stage) <- flox.stage$cluster.germ

    #print("preparing to find diff genes")
    flox.stage <- PrepSCTFindMarkers(flox.stage, assay = "SCT_GERM", verbose = TRUE)

    print("finding diff genes")
    diff.genes <- FindAllMarkers(flox.stage, assay = "SCT_GERM", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(flox.stage)[2]/length(levels(flox.stage$cluster.germ)), return.thresh = 0.05)

    filename.robj <- paste0("diff-genes-CC-flox-", i, "-", length(levels(flox.stage$cluster.germ)), ".Robj")
    save(diff.genes, file = filename.robj)

    filename.csv <- paste0("diff-genes-CC-flox-", i, "-", length(levels(flox.stage$cluster.germ)), ".csv")
    write.csv(x = diff.genes, file = filename.csv)
    
    rm(flox.stage)

}
