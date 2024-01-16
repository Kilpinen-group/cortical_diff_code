#!/bin/sh
#$ -S /bin/sh
#$ -N crangerCount
#$ -q seq_sca.q
#$ -p None
#$ -sync n
#$ -l mem_free=20G
#$ -l h_rt=604800
#$ -e log
#$ -o log
#$ -l stars=1

SAMPLEID=$(head -n $SGE_TASK_ID $1 | tail -n1)

## create execution file for this script
cd cranger/Batch1-2_220921-221021/count_211123_A00464_0418_BHNVTFDSX2
export MPSTKZ=8M
cellranger-6.1.2/cellranger count --id=${SAMPLEID} --transcriptome=refs/GRCh38/refdata-gex-GRCh38-2020-A/ --fastqs=211123_A00464_0418_BHNVTFDSX2/ --sample=${SAMPLEID} --expect-cells=16000

