#!/usr/bin/env bash

Rscript ./Scripts/BasicAAAlignment.R ./filtered_JackHMMER_sequences_all_final.fa
python3.7m Scripts/SeqCluster_v0.0/SeqCluster.py basicAln_output.aln
Rscript ./Scripts/BasicAAAlignment.R ./outputSeqs.fa
