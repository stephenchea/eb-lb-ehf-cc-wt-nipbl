#!/bin/bash

#SBATCH --job-name=cell_ranger      ## job name
#SBATCH -A schea2     ## account to charge
#SBATCH -p free  ## run on the standard partition
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=10       ## number of OpenMP threads
#SBATCH --error=slurm-%J.err ## error log file

module purge

module load openmpi/4.0.3/gcc.8.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd /dfs6/pub/schea2/20210505-aggr-bud

#cp -r /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/aggregated/20171003ab-true-flox-20180515ab-20181220abc-aggr-mapped/outs/filtered_feature_bc_matrix ./20171003ab-true-flox-20180515ab-20181220abc-aggr-mapped

module load R/4.0.4

#run the Rscript named ##-seurat.R and generate ##-seruat.Rhistory
R --no-save --quiet < 05a-seurat.R > 05a-seurat.Rhistory

