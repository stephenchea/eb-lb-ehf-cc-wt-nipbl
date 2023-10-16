#!/bin/bash

#SBATCH --job-name=count_20190722abc      ## job name
#SBATCH -A schea2     ## account to charge
#SBATCH -p free  ## run on the standard partition
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=40       ## number of OpenMP threads
#SBATCH --error=slurm-%J.err ## error log file

module purge

module load openmpi/4.0.3/gcc.8.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load cellranger/3.0.2

#run cellranger count for 20180515ab
cd /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515ab-count-force-cells

#remove presumptive wildtype libraries
rm -r AYR4-57

echo "running AYR4-57"
cellranger count --id=AYR4-57 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515b-raw --sample=AYR4_57,AYR4-57 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2750
