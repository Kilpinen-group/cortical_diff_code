#!/bin/bash -l
#$ -cwd
#$ -N demuxlet
#$ -q hugemem.q
#$ -l mem_free=20G
#$ -l h_rt=604800
#$ -e log
#$ -o log

INPUTBAM=$(head -n $SGE_TASK_ID $1 | tail -n1)
echo $INPUTBAM
DEMUXLET=/soft/bin/demuxlet
rootfolder=$(echo $INPUTBAM | sed 's/possorted_genome_bam.bam//g')
INPUTVCF=data/vcf/curatedGeno.recoded.vcf
INPUTBARCODES=$(echo ${rootfolder}filtered_feature_bc_matrix/barcodes.tsv.gz)
OUTPUT=$(echo ${rootfolder}output.demuxlet.doublet0.5)
SAMPLELIST=sample_list.txt

$DEMUXLET --sam $INPUTBAM --vcf $INPUTVCF --field GT --doublet-prior 0.5 --group-list $INPUTBARCODES --out $OUTPUT --sm-list $SAMPLELIST

 
