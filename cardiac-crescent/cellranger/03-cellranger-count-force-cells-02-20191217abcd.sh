#!/bin/bash

#SBATCH --job-name=count_20191217abcd      ## job name
#SBATCH -A alcalof_lab     ## account to charge
#SBATCH -p standard  ## run on the standard partition
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --cpus-per-task=40       ## number of OpenMP threads
#SBATCH --error=slurm-%J.err ## error log file

module purge

module load openmpi/4.0.3/gcc.8.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#run cellranger count for 20191217abcd
cd /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd

cd 20191217abcd-count-force-cells

module load cellranger/3.0.2

#for each sample, run cellranger count

cellranger count --id=FBQ3-7 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217d-raw --sample=FIN2,CaloA_FIN2 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=7014

cellranger count --id=EZC4-2 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217d-raw --sample=FIN3,CaloA_FIN3 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=5175

cellranger count --id=FBQ3-2 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217c-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20191217abcd/20191217d-raw --sample=FIN4,CaloA_FIN4 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=8890
