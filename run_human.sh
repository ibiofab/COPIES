#!/bin/bash

source venv38/bin/activate

# enable MKL for amd cpu
export MKL_DEBUG_CPU_TYPE=5

time python /home/ubuntu/COPIES/code/main.py -g GRCh38/GCF_000001405.40_GRCh38.p14_genomic.fna -t GRCh38/proteins_51_1820449.csv \
 -p NGG -o 3prime -l 20 -sl 10 --edit_dist 6 --gene_density_len 10000 -hr_l 50 \
 --protein_file GRCh38/GCF_000001405.40_GRCh38.p14_protein.faa --blast_org 'Homo sapiens' -out ../output.csv

echo "done!"
