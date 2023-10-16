library(Seurat)
library(tidyverse)
library(pracma)
library(RColorBrewer)
library(scales)
library(patchwork)

for (i in c(1, 2, 3)) {

print(i)
    
print("setting working directory")
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat/flox-integr-res-0.0180")

print("loading neighbors")
filename.robj <- paste0("flox-", i, "-neighbors.Robj")
load(filename.robj)

print("removing assays")
flox.clust[["SCT_GENOTYPE"]] <- NULL
flox.clust[["SCT_EMBRYO"]] <- NULL
flox.clust[["integrated"]] <- NULL

print("saving neighbors")
save(flox.clust, file = filename.robj)

}
