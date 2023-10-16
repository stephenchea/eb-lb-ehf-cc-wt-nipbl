#!/bin/bash

#SBATCH --job-name=seurat      ## job name
#SBATCH -A alcalof_lab     ## account to charge
#SBATCH -p free  ## run on the standard partition
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=10      ## number of OpenMP threads
#SBATCH --error=slurm-%J.err ## error log file

module load R/4.0.4

#move into working directory
cd /dfs6/pub/schea2/20210505-seurat

#run the Rscript named ##-seurat.R and generate ##-seruat.Rhistory
R --no-save --quiet < 47b-set-clusters.R > 47b-set-clusters.Rhistory
