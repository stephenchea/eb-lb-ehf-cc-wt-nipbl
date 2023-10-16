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
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120")

    print("loading umap")
    filename.robj <- paste0("fin-", i, "-projected.Robj")
    load(file = filename.robj)
    
    fin.germ$cluster.germ <- fin.germ$cluster.whole.flox
    fin.germ$cluster.germ <- factor(fin.germ$cluster.germ)
        
    print("creating new directory")
    name.direct <- paste0("fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ)))
    dir.create(name.direct)
    
    print("setting working directory")
    path.direct <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/", name.direct)
    setwd(path.direct)
    
    Idents(fin.germ) <- fin.germ$stage

    #subsetting cluster stage
    fin.stage <- subset(fin.germ, idents = "LB")
    
    rm(fin.germ)
    
    print("performing sctransform on entire cluster")
    fin.stage <- SCTransform(object = fin.stage, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.stage$label.embryo))), n_genes = round(mean(fin.stage$nFeature_RNA)), new.assay.name = "SCT_GERM", verbose = TRUE)
    
    Idents(fin.stage) <- fin.stage$cluster.germ

    #print("preparing to find diff genes")
    fin.stage <- PrepSCTFindMarkers(fin.stage, assay = "SCT_GERM", verbose = TRUE)

    print("finding diff genes")
    diff.genes <- FindAllMarkers(fin.stage, assay = "SCT_GERM", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.stage)[2]/length(levels(fin.stage$cluster.germ)), return.thresh = 0.05)

    filename.robj <- paste0("diff-genes-LB-fin-", i, "-projected-flox-", i, "-", length(levels(fin.stage$cluster.germ)), ".Robj")
    save(diff.genes, file = filename.robj)

    filename.csv <- paste0("diff-genes-LB-fin-", i, "-projected-flox-", i, "-", length(levels(fin.stage$cluster.germ)), ".csv")
    write.csv(x = diff.genes, file = filename.csv)
    
    rm(fin.stage)
    
    print(i)
    
    print("setting working directory")
    setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120")

    print("loading umap")
    filename.robj <- paste0("fin-", i, "-projected.Robj")
    load(file = filename.robj)
    
    fin.germ$cluster.germ <- fin.germ$cluster.whole.flox
    fin.germ$cluster.germ <- factor(fin.germ$cluster.germ)
        
    print("creating new directory")
    name.direct <- paste0("fin-", i, "-projected-flox-", i, "-", length(levels(fin.germ$cluster.germ)))
    dir.create(name.direct)
    
    print("setting working directory")
    path.direct <- paste0("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120/", name.direct)
    setwd(path.direct)
    
    Idents(fin.germ) <- fin.germ$stage

    #subsetting cluster stage
    fin.stage <- subset(fin.germ, idents = "CC")
    
    rm(fin.germ)
    
    print("performing sctransform on entire cluster")
    fin.stage <- SCTransform(object = fin.stage, assay = "RNA", variable.features.n = NULL, variable.features.rv.th = 1.1, return.only.var.genes = FALSE, min_cells = 2, n_cells = round(mean(table(fin.stage$label.embryo))), n_genes = round(mean(fin.stage$nFeature_RNA)), new.assay.name = "SCT_GERM", verbose = TRUE)
    
    Idents(fin.stage) <- fin.stage$cluster.germ

    #print("preparing to find diff genes")
    fin.stage <- PrepSCTFindMarkers(fin.stage, assay = "SCT_GERM", verbose = TRUE)

    print("finding diff genes")
    diff.genes <- FindAllMarkers(fin.stage, assay = "SCT_GERM", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.stage)[2]/length(levels(fin.stage$cluster.germ)), return.thresh = 0.05)

    filename.robj <- paste0("diff-genes-CC-fin-", i, "-projected-flox-", i, "-", length(levels(fin.stage$cluster.germ)), ".Robj")
    save(diff.genes, file = filename.robj)

    filename.csv <- paste0("diff-genes-CC-fin-", i, "-projected-flox-", i, "-", length(levels(fin.stage$cluster.germ)), ".csv")
    write.csv(x = diff.genes, file = filename.csv)
    
    rm(fin.stage)
    
}
