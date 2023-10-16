#load libraries
library(Seurat)
library(tidyverse)

#set working directory
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

#read cell ranger data
all.aggr <- Read10X(data.dir = "20171003ab-true-flox-20180515ab-20181220abc-20190722abc-20191217abcd-aggr-mapped")

#create vector for stage
stage <- vector()
stage[1:37792] <- "LB"
stage[37793:109784] <- "CC"

#factor stage
stage <- factor(x = stage, levels = c("LB", "CC"))

#create vector for collection date
date.collection <- vector()
date.collection[1:2596] <- "20171003"
date.collection[2597:6749] <- "20180515"
date.collection[6750:16973] <- "20181220"
date.collection[16974:20511] <- "20171003"
date.collection[20512:29061] <- "20180515"
date.collection[29062:37792] <- "20181220"
date.collection[37793:48473] <- "20190722"
date.collection[48474:82413] <- "20191217"
date.collection[82414:91674] <- "20190722"
date.collection[91675:109784] <- "20191217"

#factor collection
date.collection <- factor(x = date.collection, levels = c("20171003", "20180515", "20181220", "20190722", "20191217"))

#create vector for collection label
label.collection <- vector()
label.collection[1:2596] <- "CL1"
label.collection[2597:6749] <- "CL2"
label.collection[6750:16973] <- "CL3"
label.collection[16974:20511] <- "CL1"
label.collection[20512:29061] <- "CL2"
label.collection[29062:37792] <- "CL3"
label.collection[37793:48473] <- "CL4"
label.collection[48474:82413] <- "CL5"
label.collection[82414:91674] <- "CL4"
label.collection[91675:109784] <- "CL5"

label.collection <- factor(x = label.collection, levels = c("CL1", "CL2", "CL3", "CL4", "CL5"))

#create vector for embryo id
id.embryo <- vector()
id.embryo[1:2596] <- "AVJ4-56"
id.embryo[2597:3999] <- "AYQ8-11"
id.embryo[4000:6749] <- "AYR4-57"
id.embryo[6750:14351] <- "BBH6-17"
id.embryo[14352:16973] <- "AXL7-21"
id.embryo[16974:18958] <- "ATU6-45"
id.embryo[18959:20511] <- "ATU6-69"
id.embryo[20512:24972] <- "AYR4-22"
id.embryo[24973:29061] <- "AZK2-45"
id.embryo[29062:32724] <- "BBH6-59"
id.embryo[32725:37792] <- "AXL7-69"
id.embryo[37793:39912] <- "EQQ5-4"
id.embryo[39913:43676] <- "ERU4-4"
id.embryo[43677:48473] <- "ERU4-7"
id.embryo[48474:54146] <- "EYO8-1"
id.embryo[54147:61968] <- "EYO8-6"
id.embryo[61969:70306] <- "EZC4-8"
id.embryo[70307:75399] <- "FBQ3-4"
id.embryo[75400:82413] <- "FBQ3-7"
id.embryo[82414:84705] <- "EQQ5-2"
id.embryo[84706:87005] <- "ERZ4-4"
id.embryo[87006:87734] <- "ERZ4-6"
id.embryo[87735:89479] <- "ERZ4-7"
id.embryo[89480:91674] <- "ESH7-1"
id.embryo[91675:95719] <- "EYO8-3"
id.embryo[95720:100894] <- "EZC4-2"
id.embryo[100895:109784] <- "FBQ3-2"

#factor embryo id
id.embryo <- factor(x = id.embryo, levels = c("AVJ4-56",
"AYQ8-11",
"AYR4-57",
"BBH6-17",
"AXL7-21",
"ATU6-45",
"ATU6-69",
"AYR4-22",
"AZK2-45",
"BBH6-59",
"AXL7-69",
"EQQ5-4",
"ERU4-4",
"ERU4-7",
"EYO8-1",
"EYO8-6",
"EZC4-8",
"FBQ3-4",
"FBQ3-7",
"EQQ5-2",
"ERZ4-4",
"ERZ4-6",
"ERZ4-7",
"ESH7-1",
"EYO8-3",
"EZC4-2",
"FBQ3-2"))

#create vector for embryo label
label.embryo <- vector()
label.embryo[1:2596] <- "LB1"
label.embryo[2597:3999] <- "LB2"
label.embryo[4000:6749] <- "LB3"
label.embryo[6750:14351] <- "LB4"
label.embryo[14352:16973] <- "LB5"
label.embryo[16974:18958] <- "LB6"
label.embryo[18959:20511] <- "LB7"
label.embryo[20512:24972] <- "LB8"
label.embryo[24973:29061] <- "LB9"
label.embryo[29062:32724] <- "LB10"
label.embryo[32725:37792] <- "LB11"
label.embryo[37793:39912] <- "CC1"
label.embryo[39913:43676] <- "CC2"
label.embryo[43677:48473] <- "CC3"
label.embryo[48474:54146] <- "CC4"
label.embryo[54147:61968] <- "CC5"
label.embryo[61969:70306] <- "CC6"
label.embryo[70307:75399] <- "CC7"
label.embryo[75400:82413] <- "CC8"
label.embryo[82414:84705] <- "CC9"
label.embryo[84706:87005] <- "CC10"
label.embryo[87006:87734] <- "CC11"
label.embryo[87735:89479] <- "CC12"
label.embryo[89480:91674] <- "CC13"
label.embryo[91675:95719] <- "CC14"
label.embryo[95720:100894] <- "CC15"
label.embryo[100895:109784] <- "CC16"

#factor embryo label
label.embryo <- factor(x = label.embryo, levels = c("LB1",
"LB2",
"LB3",
"LB4",
"LB5",
"LB6",
"LB7",
"LB8",
"LB9",
"LB10",
"LB11",
"CC1",
"CC2",
"CC3",
"CC4",
"CC5",
"CC6",
"CC7",
"CC8",
"CC9",
"CC10",
"CC11",
"CC12",
"CC13",
"CC14",
"CC15",
"CC16"))

#create vector for genotype
genotype <- vector()
genotype[1:16973] <- "Nipbl Flox/+"
genotype[16974:37792] <- "Nipbl FIN/+"
genotype[37793:82413] <- "Nipbl Flox/+"
genotype[82414:109784] <- "Nipbl FIN/+"

#factor genotype
genotype <- factor(x = genotype, levels = c("Nipbl Flox/+", "Nipbl FIN/+"))

#create metadata
metadata <- data.frame(stage, date.collection, label.collection, id.embryo, label.embryo, genotype)
rownames(metadata) <- colnames(all.aggr)

#create seurat object
all.aggr <- CreateSeuratObject(counts = all.aggr, project = "all.aggr", min.cells = 0, min.features = 0, meta.data = metadata)

#split object by embryo
embryo.list <- SplitObject(object = all.aggr, split.by = "label.embryo")

#create pseudobulk
umi.embryo <- lapply(X = embryo.list, FUN = function(x) rowSums(x@assays$RNA@counts))

#aggregate pseudobulk
umi.embryo <- data.frame(umi.embryo[[1]], umi.embryo[[2]], umi.embryo[[3]], umi.embryo[[4]], umi.embryo[[5]], umi.embryo[[6]], umi.embryo[[7]], umi.embryo[[8]], umi.embryo[[9]], umi.embryo[[10]], umi.embryo[[11]], umi.embryo[[12]], umi.embryo[[13]], umi.embryo[[14]], umi.embryo[[15]], umi.embryo[[16]], umi.embryo[[17]], umi.embryo[[18]], umi.embryo[[19]], umi.embryo[[20]], umi.embryo[[21]], umi.embryo[[22]], umi.embryo[[23]], umi.embryo[[24]], umi.embryo[[25]], umi.embryo[[26]], umi.embryo[[27]])

#find genes that are expressed in at least one embryo
detected.genes <- apply(X = umi.embryo, MARGIN = 1, FUN = function(x) sum(x) > 0)

#keep only those genes from this point foward
all.aggr <- subset(x = all.aggr, features = rownames(all.aggr)[detected.genes])

save(all.aggr, file = "all-aggr-seurat.Robj")
