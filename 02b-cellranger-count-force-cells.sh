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

#run cellranger count for 20181220abc
cd /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220abc-count-force-cells

#remove presumptive wildtype libraries
rm -r BBH6-17
rm -r AXL7-21

echo "running BBH6-17"
#for each sample, run cellranger count
cellranger count --id=BBH6-17 --fastqs=/dfs3b/alcalof_lab/lander-calof/incoming/20181220abc/20181220a-raw,/dfs3b/alcalof_lab/lander-calof/incoming/20181220abc/20181220b-raw,/dfs3b/alcalof_lab/lander-calof/incoming/20181220abc/20181220c-raw --sample=W1 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=7602

echo "running AXL7-21"
cellranger count --id=AXL7-21 --fastqs=/dfs3b/alcalof_lab/lander-calof/incoming/20181220abc/20181220a-raw,/dfs3b/alcalof_lab/lander-calof/incoming/20181220abc/20181220b-raw,/dfs3b/alcalof_lab/lander-calof/incoming/20181220abc/20181220c-raw --sample=W2 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2622
