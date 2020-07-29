#!/usr/bin/env bash

inputFile=$1

echo $inputFile
K=$2

Method="BIRCH"
BASR_PATH="/Users/ptdolan/Documents/GitHub/BASR"

Rscript ${BASR_PATH}/Scripts/BasicAAAlignment.R $inputFile
python3.7m ${BASR_PATH}/Scripts/SeqCluster_v0.0/SeqCluster.py ${inputFile}_basicAA.aln $K $Method
Rscript ${BASR_PATH}/Scripts/BasicAAAlignment.R ${inputFile}_basicAA.aln${Method}_outputSeqs.fasta
cat ${inputFile}_basicAA.aln${Method}_outputSeqs.fasta_basicAA.aln > ${inputFile/.fa*/}_${K}reps.aln
sed -i'' -e 's/(/_/g' ${inputFile/.fa*/}_${K}reps.aln
sed -i'' -e 's/)/_/g' ${inputFile/.fa*/}_${K}reps.aln
sed -i'' -e 's/;/_/g' ${inputFile/.fa*/}_${K}reps.aln
sed -i'' -e 's/:/_/g' ${inputFile/.fa*/}_${K}reps.aln
sed -i'' -e 's/,/_/g' ${inputFile/.fa*/}_${K}reps.aln
FastTree -gtr ${inputFile/.fa*/}_${K}reps.aln > ${inputFile/.fa*/}_${K}reps.tree
echo Alignment reduced to $K representatives.
