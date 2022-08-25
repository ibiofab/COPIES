#!/bin/bash

source venv38/bin/activate

time python code/main.py -g s288c/GCF_000146045.2_R64_genomic.fna -t s288c/proteins_15_22535.csv \
 -p NGG -o 3prime -l 20 -sl 10 --edit_dist 6 --gene_density_len 10000 -hr_l 50 \
 --protein_file s288c/GCF_000146045.2_R64_protein.faa --blast_org 'Saccharomyces cerevisiae' -out ../output.csv

echo "done!"
