#!/bin/bash

#SBATCH --job-name=count_20190722abc      ## job name
#SBATCH -A alcalof_lab     ## account to charge
#SBATCH -p free  ## run on the standard partition
#SBATCH --nodes=1               ## number of nodes the job will use
#SBATCH --ntasks=1              ## number of processes to launch
#SBATCH --mem=240G              ## request 100GB of memory
#SBATCH --cpus-per-task=40       ## number of OpenMP threads
#SBATCH --error=slurm-%J.err ## error log file

module purge

module load openmpi/4.0.3/gcc.8.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#run cellranger count for 20191217abcd
cd /share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722abc-count-force-cells

module load cellranger/3.0.2

#for each sample, run cellranger count
#cellranger count --id=ERU4-4 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=FLOX_CC1,FLOW_CC1 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=3764

#cellranger count --id=ERU4-7 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=FLOX_CC2,FLOW_CC2 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=4797

#cellranger count --id=EQQ5-2 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=FLOX_CC3,FLOW_CC3 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2292

cellranger count --id=EQQ5-4 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=FLOX_CC4,FLOW_CC4 --transcriptome=/dfs6/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2120

#cellranger count --id=ERZ4-4 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=Fin_CC1,FIN_CC1 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2300

#cellranger count --id=ERZ4-6 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw --sample=Fin_CC2,FIN_CC2 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=729

#cellranger count --id=ERZ4-7 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=Fin_CC3,FIN_CC3 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=1745

#cellranger count --id=ESH7-1 --fastqs=/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722a-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722b-raw,/share/crsp/lab/alcalof/share/dfs3-backup/lander-calof/incoming/20190722abc/20190722c-raw --sample=Fin_CC4,FIN_CC4 --transcriptome=/pub/schea2/refs/mod-mm10/mod-mm10 --force-cells=2195
