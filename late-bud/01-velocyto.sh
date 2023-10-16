#!/bin/bash
#$ -N velocyto
#$ -q bio,math,bme,pub64,abio,free64
#$ -pe openmp 8-64
#$ -ckpt restart
#$ -m beas

#move into general directory containing 10X outputs
cd /pub/schea2/20200824

module load anaconda
PATH=/data/users/schea2/.local/bin:$PATH

#generate looms for 20171003ab

echo "generating looms for 20171003ab"

for i in "ATU6-45" "ATU6-69" "AVJ4-23" "AVJ4-56"

do

echo "copying input files for"
echo $i

cd /dfs5/bio/lander-calof/incoming/20171003ab/20171003ab-count-force-cells/$i/outs/filtered_feature_bc_matrix

cp barcodes.tsv.gz /pub/schea2/20200824

cd /dfs5/bio/lander-calof/incoming/20171003ab/20171003ab-count-force-cells/$i/outs

cp possorted_genome_bam.bam /pub/schea2/20200824

cd /pub/schea2/20200824

gunzip barcodes.tsv.gz

echo "generating loom for"
echo $i

velocyto run -b barcodes.tsv -o /pub/schea2/20200824 -e $i -m mm10_rmsk.gtf possorted_genome_bam.bam /pub/schea2/refs/mod-mm10/mod-mm10-genes.gtf

rm barcodes.tsv
rm possorted_genome_bam.bam
rm cellsorted_possorted_genome_bam.bam

echo "removing input files for"
echo $i
done

#generate looms for 20180515ab

echo "generating looms for 20180515ab"

for i in "AYQ8-11" "AYR4-22" "AYR4-57" "AZK2-45"

do

echo "copying input files for"
echo $i

cd /dfs5/bio/lander-calof/incoming/20180515ab/20180515ab-count-force-cells/$i/outs/filtered_feature_bc_matrix

cp barcodes.tsv.gz /pub/schea2/20200824

cd /dfs5/bio/lander-calof/incoming/20180515ab/20180515ab-count-force-cells/$i/outs

cp possorted_genome_bam.bam /pub/schea2/20200824

cd /pub/schea2/20200824

gunzip barcodes.tsv.gz

echo "generating loom for"
echo $i

velocyto run -b barcodes.tsv -o /pub/schea2/20200824 -e $i -m mm10_rmsk.gtf possorted_genome_bam.bam /pub/schea2/refs/mod-mm10/mod-mm10-genes.gtf

rm barcodes.tsv
rm possorted_genome_bam.bam
rm cellsorted_possorted_genome_bam.bam

echo "removing input files for"
echo $i
done

#generate looms for 20181220abc

echo "generating looms for 20181220abc"

for i in "AXL7-21" "AXL7-69" "BBH6-17" "BBH6-59"

do

echo "copying input files for"
echo $i

cd /dfs5/bio/lander-calof/incoming/20181220abc/20181220abc-count-force-cells/$i/outs/filtered_feature_bc_matrix

cp barcodes.tsv.gz /pub/schea2/20200824

cd /dfs5/bio/lander-calof/incoming/20181220abc/20181220abc-count-force-cells/$i/outs

cp possorted_genome_bam.bam /pub/schea2/20200824

cd /pub/schea2/20200824

gunzip barcodes.tsv.gz

echo "generating loom for"
echo $i

velocyto run -b barcodes.tsv -o /pub/schea2/20200824 -e $i -m mm10_rmsk.gtf possorted_genome_bam.bam /pub/schea2/refs/mod-mm10/mod-mm10-genes.gtf

rm barcodes.tsv
rm possorted_genome_bam.bam
rm cellsorted_possorted_genome_bam.bam

echo "removing input files for"
echo $i
done
