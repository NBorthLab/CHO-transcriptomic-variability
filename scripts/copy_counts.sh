#!/usr/bin/env bash

: '
Copy the counts matrix from the HPC cluster to the local machine.
'

rsync -avzhPR rkn07:/data/borth/mriedl/projects/P06-GeneVariability/./results/counts/all_gene_counts.rds .
