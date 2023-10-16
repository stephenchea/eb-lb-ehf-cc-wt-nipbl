#!/bin/bash

#SBATCH --job-name=aggr_20190722abc_20191217abcd      ## job name
#SBATCH -A alcalof_lab     ## account to charge
#SBATCH -p standard  ## run on the standard partition
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=30       ## number of OpenMP threads
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --mail-type=fail,end
#SBATCH --mail-user=schea2@uci.edu

module purge

module load openmpi/4.0.3/gcc.8.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#move into the folder in which we want to store our aggregated matrix
cd /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/aggregated

#load cell ranger
module load cellranger/3.0.2

#run cellranger aggr on samples
cellranger aggr --id=20190722abc-20191217abcd-aggr-mapped --csv=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/aggregated/20190722abc-20191217abcd-libs.csv --normalize=mapped
