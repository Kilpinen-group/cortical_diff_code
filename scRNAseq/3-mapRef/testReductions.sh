#!/bin/bash
  
echo $LSB_JOBINDEX
echo $1
line=$(head -n $LSB_JOBINDEX $1 | tail -n1)
queryRed=$(echo $line | cut -f1 -d" ")
refRed=$(echo $line | cut -f2 -d" ")
dataset=$(echo $line | cut -f3 -d" ")

echo "Starting mapping reference"
echo Q${queryRed} R${refRed} D${dataset}
Rscript testReductions.R $queryRed $refRed $dataset

