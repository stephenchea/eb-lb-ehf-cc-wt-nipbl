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

#run cellranger count for 20171003ab
/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20171003ab/20171003ab-count-force-cells

#remove presumptive wildtype libraries
rm -r AVJ4-23
rm -r AVJ4-56

echo "running AVJ4-23"
#for each sample, run cellranger count
cellranger count --id=AVJ4-23 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/heart-rnaseq-10x-oct-2017/20171003a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/heart-rnaseq-10x-oct-2017/20171003b-raw --sample=AvJ4-23,AVJ4-23 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=3945

echo "running AVJ4-56"
cellranger count --id=AVJ4-56 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/heart-rnaseq-10x-oct-2017/20171003a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/heart-rnaseq-10x-oct-2017/20171003b-raw --sample=AvJ4-56,AVJ4-56 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2596

#run cellranger count for 20180515ab
/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515ab-count-force-cells

#remove presumptive wildtype libraries
rm -r AYQ8-11
rm -r AYR4-57

echo "running AYQ8-11"
#for each sample, run cellranger count
cellranger count --id=AYQ8-11 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515b-raw --sample=AYQ8_11,AYQB-11 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=1403

echo "running AYR4-57"
cellranger count --id=AYR4-57 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20180515ab/20180515b-raw --sample=AYR4_57 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2750

#run cellranger count for 20181220abc
/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220abc-count-force-cells

#remove presumptive wildtype libraries
rm -r BBH6-17
rm -r AXL7-21

echo "running BBH6-17"
#for each sample, run cellranger count
cellranger count --id=BBH6-17 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220c-raw --sample=W1 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=7602

echo "running AXL7-21"
cellranger count --id=AXL7-21 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20181220abc/20181220c-raw --sample=W2 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2622
