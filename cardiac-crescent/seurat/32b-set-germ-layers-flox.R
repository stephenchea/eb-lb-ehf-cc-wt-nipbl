library(Seurat)
library(tidyverse)

#set working directory
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

#load seurat object
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-neighbors.Robj")

#set working directory
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120")

#load cluster identities
load("umap-by-cluster-flox-integr-res-0.0120-26.Robj")
flox.integr$cluster.whole <- umap.whole$cluster.whole

#load germ layer identities
load("umap-by-germ-layer-flox-integr-res-0.0120-26.Robj")
flox.integr$germ.whole <- umap.whole$germ.whole

#set identity to germ layer
Idents(flox.integr) <- umap.whole$germ.whole

#print("preparing to find diff genes")
#flox.integr <- PrepSCTFindMarkers(flox.integr, assay = "SCT_GENOTYPE", verbose = TRUE)

#print("finding diff genes")
#diff.genes <- FindAllMarkers(flox.integr, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(flox.integr)[2]/length(levels(flox.integr$germ.whole)), return.thresh = 0.05)

#save(diff.genes, file = "diff-genes-between-germ-layers-flox-integr-res-0.0120-26.Robj")

#write.csv(x = diff.genes, file = "diff-genes-between-germ-layers-flox-integr-res-0.0120-26.csv")

############################################

#set identity to germ layer
Idents(flox.integr) <- flox.integr$stage

#subsetting cluster stage
flox.stage <- subset(flox.integr, idents = "LB")

rm(flox.integr)

#set identity to germ layer
Idents(flox.stage) <- flox.stage$germ.whole

#print("preparing to find diff genes")
#flox.stage <- PrepSCTFindMarkers(flox.stage, assay = "SCT_GENOTYPE", verbose = TRUE)

print("finding diff genes")
diff.genes <- FindAllMarkers(flox.stage, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(flox.stage)[2]/length(levels(flox.stage$germ.whole)), return.thresh = 0.05)

save(diff.genes, file = "diff-genes-between-germ-layers-LB-flox-integr-res-0.0120-26.Robj")

write.csv(x = diff.genes, file = "diff-genes-between-germ-layers-LB-flox-integr-res-0.0120-26.csv")

rm(flox.stage)

############################################

#set working directory
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

#load seurat object
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-neighbors.Robj")

#set working directory
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0120")

#load cluster identities
load("umap-by-cluster-flox-integr-res-0.0120-26.Robj")
flox.integr$cluster.whole <- umap.whole$cluster.whole

#load germ layer identities
load("umap-by-germ-layer-flox-integr-res-0.0120-26.Robj")
flox.integr$germ.whole <- umap.whole$germ.whole

#set identity to germ layer
Idents(flox.integr) <- flox.integr$stage

#subsetting cluster stage
flox.stage <- subset(flox.integr, idents = "CC")

#set identity to germ layer
Idents(flox.stage) <- flox.stage$germ.whole

#print("preparing to find diff genes")
#flox.stage <- PrepSCTFindMarkers(flox.stage, assay = "SCT_GENOTYPE", verbose = TRUE)

print("finding diff genes")
diff.genes <- FindAllMarkers(flox.stage, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(flox.stage)[2]/length(levels(flox.stage$germ.whole)), return.thresh = 0.05)

save(diff.genes, file = "diff-genes-between-germ-layers-CC-flox-integr-res-0.0120-26.Robj")

write.csv(x = diff.genes, file = "diff-genes-between-germ-layers-CC-flox-integr-res-0.0120-26.csv")

