#a new folder was created at the following path
#move into this folder; it will house further analysis
setwd("/dfs6/pub/schea2/20210505-aggr-bud")

library(Seurat)
library(tidyverse)

print("loaded R libraries")

#the working directory has a folder called 20171003ab-true-flox-20180515ab-20181220abc-aggr-mapped that contains cell ids, genes, and read counts
#load cellranger aggregated dataset
all.aggr <- Read10X(data.dir = "20171003ab-true-flox-20180515ab-20181220abc-aggr-mapped")
print("assembled 10X library")

print("dimensions")
dim(all.aggr)
#record number of cells and genes in matrix

experiment <- vector()
experiment[1:2596] <- "20171003"
experiment[2597:6749] <- "20180515"
experiment[6750:16973] <- "20181220"
experiment[16974:20511] <- "20171003"
experiment[20512:29061] <- "20180515"
experiment[29062:37792] <- "20181220"

embryo.id <- vector()
embryo.id[1:2596] <- "AVJ4-56"
embryo.id[2597:3999] <- "AYQ8-11"
embryo.id[4000:6749] <- "AYR4-57"
embryo.id[6750:14351] <- "BBH6-17"
embryo.id[14352:16973] <- "AXL7-21"
embryo.id[16974:18958] <- "ATU6-45"
embryo.id[18959:20511] <- "ATU6-69"
embryo.id[20512:24972] <- "AYR4-22"
embryo.id[24973:29061] <- "AZK2-45"
embryo.id[29062:32724] <- "BBH6-59"
embryo.id[32725:37792] <- "AXL7-69"

genotype <- vector()
genotype[1:16973] <- "Flox"
genotype[16974:37792] <- "FIN"

#generate metadata of aggr.raw for aggr.seurat indicating experiment, sample.id, and genotype of each cell
all.aggr.metadata <- data.frame(experiment, embryo.id, genotype)
rownames(all.aggr.metadata) <- colnames(all.aggr)

print("generated metadata")

#Create seurat object; only keep genes that are expressed in minimum of 3 cells
all.aggr <- CreateSeuratObject(counts = all.aggr, project = "all.aggr", min.cells = 3, min.features = 0, meta.data = all.aggr.metadata)
print("generated seurat object")

print("dimensions")
dim(all.aggr)

save(all.aggr, file = "all-aggr-seurat.Robj")
