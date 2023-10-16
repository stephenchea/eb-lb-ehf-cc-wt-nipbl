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

#move into the folder in which we want to store our aggregated matrix
cd /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/aggregated

#load cell ranger
module load cellranger/3.0.2

#run cellranger aggr on samples
cellranger aggr --id=20190211ab-ehf-aggr-mapped --csv=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/aggregated/20190211ab-ehf-libs.csv --normalize=mapped
