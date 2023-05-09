# CRISPR-COPIES
**CO**mputational **P**ipeline for the **I**dentification of CRISPR/Cas-facilitated int**E**gration **S**ites (`CRISPR-COPIES`) is a command line and web interface tool for rapid discovery of neutral integration sites. Designed to work for any organism with a genome in NCBI and for any CRISPR system, `CRISPR-COPIES` can identify neutral sites in a genome-wide manner. The identified sites can be used for characterization of synthetic biology toolkits, rapid strain construction to produce valuable biochemicals, and human gene and cell therapy.

- [Installation](#installation)
- [Usage](#usage)
- [Parameters](#parameters)

### READ THE PAPER!

This repository accompanies the work ["CRISPR-COPIES: An in silico platform for discovery of neutral integration sites for CRISPR/Cas-facilitated gene integration"](https://www.google.com).

<details>
<summary>If you use this tool, please cite us:</summary>

```bibtex

```
</details>

### Installation
```
pip install -r requirements38.txt
```
### Usage

#### Command Line
`CRISPR-COPIES` can be accessed from the command line. For information on optional arguments and flags, run
```
python main.py -h
```
The parameters are explained in more details in the 'Parameter Dictionary' below. A sample example to run the script - 
```
python code/main.py -g s288c/GCF_000146045.2_R64_genomic.fna -t s288c/GCF_000146045.2_R64_feature_table.txt -p NGG -o 3prime -l 20 -sl 10 --edit_dist 6 --intspace 350 -out s288c/output.csv --distal_end_len 10000 -hr_l 50 --protein_file s288c/GCF_000146045.2_R64_protein.faa
```

#### Web Interface
You can also use `CRISPR-COPIES` through our web interface. Visit us at [CRISPR-COPIES](https://biofoundry.web.illinois.edu/copies/). 

Note: We have restricted the web interface to prokaryotic and small eukaryotic genomes. We advise you to use our command line option for genomes greater than 120 Mb in size as significant time and computation resources are required. 

### Parameters
#### The parameter dictionary for `CRISPR-COPIES` is divided into 4 sections: 
#### 1. Guide RNA
<img src=https://user-images.githubusercontent.com/60017121/175441737-28ad07ab-8888-49ec-bec3-a77c4f153292.png width="648">

#### 2. Homology Arm
<img src=https://user-images.githubusercontent.com/60017121/172052679-150a321b-be90-4d4e-939c-5233a0775ea3.png width="648">

#### 3. Harbor Information 
<img src=https://user-images.githubusercontent.com/60017121/172052744-d0394ec2-b84e-498b-a583-ca04c919b530.png width="648">

#### 4. Essentiality Information. 
<img src=https://user-images.githubusercontent.com/60017121/172052754-3719f39a-021e-42ec-b828-ec22ed6ee6a6.png width="648">
