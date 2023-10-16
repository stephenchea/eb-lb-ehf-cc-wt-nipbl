#load libraries
library(Seurat)
library(tidyverse)
library(ggsignif)
library(Rmisc)
library(scales)

#set working directory
setwd("/share/crsp/lab/alcalof/schea2/20220414-seurat")

#load cells
load("/share/crsp/lab/alcalof/schea2/20220414-seurat/all-aggr-neighbors.Robj")

#set identity to embryo
Idents(all.aggr) <- all.aggr$label.embryo

#split cells by embryo
embryo.list <- SplitObject(all.aggr, split.by = "ident")

#grab normalized counts
umi.embryo <- lapply(X = embryo.list, FUN = function(x) rowSums(x@assays$SCT_ALL@counts))

#organize counts as data frame
umi.data <- data.frame(umi.embryo[[1]], umi.embryo[[2]], umi.embryo[[3]], umi.embryo[[4]], umi.embryo[[5]], umi.embryo[[6]], umi.embryo[[7]], umi.embryo[[8]], umi.embryo[[9]], umi.embryo[[10]], umi.embryo[[11]], umi.embryo[[12]], umi.embryo[[13]], umi.embryo[[14]], umi.embryo[[15]], umi.embryo[[16]], umi.embryo[[17]], umi.embryo[[18]], umi.embryo[[19]], umi.embryo[[20]], umi.embryo[[21]], umi.embryo[[22]], umi.embryo[[23]], umi.embryo[[24]], umi.embryo[[25]], umi.embryo[[26]], umi.embryo[[27]])

#set column names to embryo
colnames(umi.data) <- names(umi.embryo)

#calculate normalization factor
cell.embryo <- table(all.aggr$label.embryo)
norm.factor <- cell.embryo / min(cell.embryo)

#divide each gene by normalization factor
norm.umi <- as.data.frame(apply(X = umi.data, MARGIN = 1, FUN = function(x) x/norm.factor))

#add stage into data frame
norm.umi$stage <- c("LB", "LB", "LB", "LB", "LB", "LB", "LB", "LB", "LB", "LB", "LB", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC")
norm.umi$stage <- factor(norm.umi$stage, levels = c("LB", "CC"))

#add genotype into data frame
norm.umi$genotype <- c("Flox", "Flox", "Flox", "Flox", "Flox", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "Flox", "Flox", "Flox", "Flox", "Flox", "Flox", "Flox", "Flox", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN", "FIN")
norm.umi$genotype <- factor(norm.umi$genotype, levels = c("Flox", "FIN"))

save(norm.umi, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/norm-umi-per-embryo.Robj")

#calculate mean and confidence interval
error.nipbl <- summarySE(data = norm.umi, measurevar = "Nipbl", groupvars = c("stage", "genotype"))

#perform t test
t.test(formula = log10(norm.umi$Nipbl) ~ norm.umi$genotype, alternative = "two.sided", paired = FALSE, conf.level = 0.95)

save(error.nipbl, file = "/share/crsp/lab/alcalof/schea2/20220414-seurat/error-nipbl-per-genotype.Robj")
