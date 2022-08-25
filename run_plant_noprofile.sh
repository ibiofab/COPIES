#!/bin/bash

source venv38/bin/activate

# enable MKL for amd cpu
export MKL_DEBUG_CPU_TYPE=5

time python /home/ubuntu/COPIES/code/main.py  -g tair10/GCF_000001735.4_TAIR10.1_genomic.fna -t tair10/GCF_000001735.4_TAIR10.1_feature_table.txt \
 -p NGG -o 3prime -l 20 -sl 10 --edit_dist 6 --gene_density_len 10000 -hr_l 50 \
 --protein_file tair10/GCF_000001735.4_TAIR10.1_protein.faa --blast_org 'Arabidopsis thaliana' -out ../output.csv --on_target 'sgRNA_ecoli(Cas9)' 

echo "done!"
