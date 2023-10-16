library(Seurat)
library(tidyverse)

# print("setting working directory")
# setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120")
# 
# print("loading fin cells")
# load("fin-integr-projected-flox-integr-res-0.0120.Robj")
# 
# print("setting fin projected cluster identities")
# load("projected-umap-by-stage-embryo-cluster-fin-integr-projected-flox-integr-res-0.0120-26.Robj")
# fin.integr$cluster.whole.flox <- whole.umap$cluster.whole.flox
# 
# print("setting fin projected germ layer identities")
# load("projected-umap-by-germ-layer-fin-integr-projected-flox-integr-res-0.0120-26.Robj")
# fin.integr$germ.whole.flox <- whole.umap$germ.whole.flox
# 
# print("setting germ layer identity")
# Idents(fin.integr) <- fin.integr$germ.whole.flox
# 
# #print("preparing to find diff genes")
# #fin.integr <- PrepSCTFindMarkers(fin.integr, assay = "SCT_GENOTYPE", verbose = TRUE)
# 
# print("finding diff genes")
# diff.genes <- FindAllMarkers(fin.integr, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.integr)[2]/length(levels(fin.integr$germ.whole.flox)), return.thresh = 0.05)
# 
# save(diff.genes, file = "diff-genes-between-germ-layers-fin-integr-projected-flox-integr-res-0.0120-26.Robj")
# 
# write.csv(x = diff.genes, file = "diff-genes-between-germ-layers-fin-integr-projected-flox-integr-res-0.0120-26.csv")
# 
# ############################################
# 
# #set identity to germ layer
# Idents(fin.integr) <- fin.integr$stage
# 
# #subsetting cluster stage
# fin.stage <- subset(fin.integr, idents = "LB")
# 
# rm(fin.integr)
# 
# #set identity to germ layer
# Idents(fin.stage) <- fin.stage$germ.whole.flox
# 
# #print("preparing to find diff genes")
# #fin.stage <- PrepSCTFindMarkers(fin.stage, assay = "SCT_GENOTYPE", verbose = TRUE)
# 
# print("finding diff genes")
# diff.genes <- FindAllMarkers(fin.stage, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.stage)[2]/length(levels(fin.stage$germ.whole.flox)), return.thresh = 0.05)
# 
# save(diff.genes, file = "diff-genes-between-germ-layers-LB-fin-integr-projected-flox-integr-res-0.0120-26.Robj")
# 
# write.csv(x = diff.genes, file = "diff-genes-between-germ-layers-LB-fin-integr-projected-flox-integr-res-0.0120-26.csv")
# 
# rm(fin.stage)

############################################

print("setting working directory")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/fin-integr-projected-flox-integr-res-0.0120")

print("loading fin cells")
load("fin-integr-projected-flox-integr-res-0.0120.Robj")

print("setting fin projected cluster identities")
load("projected-umap-by-stage-embryo-cluster-fin-integr-projected-flox-integr-res-0.0120-26.Robj")
fin.integr$cluster.whole.flox <- whole.umap$cluster.whole.flox

print("setting fin projected germ layer identities")
load("projected-umap-by-germ-layer-fin-integr-projected-flox-integr-res-0.0120-26.Robj")
fin.integr$germ.whole.flox <- whole.umap$germ.whole.flox

#set identity to germ layer
Idents(fin.integr) <- fin.integr$stage

#subsetting cluster stage
fin.stage <- subset(fin.integr, idents = "CC")

rm(fin.integr)

#set identity to germ layer
Idents(fin.stage) <- fin.stage$germ.whole.flox

#print("preparing to find diff genes")
#fin.stage <- PrepSCTFindMarkers(fin.stage, assay = "SCT_GENOTYPE", verbose = TRUE)

print("finding diff genes")
diff.genes <- FindAllMarkers(fin.stage, assay = "SCT_GENOTYPE", logfc.threshold = 0, test.use = "wilcox", slot = "data", min.pct = 0, verbose = TRUE, only.pos = TRUE, max.cells.per.ident = dim(fin.stage)[2]/length(levels(fin.stage$germ.whole.flox)), return.thresh = 0.05)

save(diff.genes, file = "diff-genes-between-germ-layers-CC-fin-integr-projected-flox-integr-res-0.0120-26.Robj")

write.csv(x = diff.genes, file = "diff-genes-between-germ-layers-CC-fin-integr-projected-flox-integr-res-0.0120-26.csv")
