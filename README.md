# COPIES
**CO**mputational **P**ipeline for **I**dentification of CRISPR-facilitated int**E**gration **S**ites (`COPIES`) is a command line and web interface tool for rapid discovery of genomic safe harbors. Designed to work for any organism with a genome in NCBI and for any CRISPR system, `COPIES` can identify neutral sites in a genome-wide manner. The identified sites can be used for synthetic biology toolkit characterization and construction of genetically stable strains for biochemical production.

- [Installation](#installation)
- [Usage](#usage)
- [Parameters](#parameters)

### READ THE PAPER!

This repository accompanies the work ["COPIES, an in-silico platform for screening safe harbors for CRISPR-facilitated gene integration"](https://www.google.com).

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
`COPIES` can be accessed from the command line. For information on optional arguments and flags, run
```
python sites.py -h
```
The parameters are explained in more details in the 'Parameter Dictionary' below. A sample example to run the script - 
```
python sites.py -
```

#### Web Interface
You can also use `COPIES` through our web interface. Visit us at [COPIES](https://www.google.com). 

Note: We have restricted the web interface to prokaryotic and small eukaryotic genomes. We advise you to use our command line option for genomes greater than __ in size as significant time and computation resources are required. 

### Parameters
#### The parameter dictionary for `COPIES` is divided into 4 sections: 
#### 1. Guide RNA
<img src=https://user-images.githubusercontent.com/60017121/175441737-28ad07ab-8888-49ec-bec3-a77c4f153292.png width="648">

#### 2. Homology Arm
<img src=https://user-images.githubusercontent.com/60017121/172052679-150a321b-be90-4d4e-939c-5233a0775ea3.png width="648">

#### 3. Harbor Information 
<img src=https://user-images.githubusercontent.com/60017121/172052744-d0394ec2-b84e-498b-a583-ca04c919b530.png width="648">

#### 4. Essentiality Information. 
<img src=https://user-images.githubusercontent.com/60017121/172052754-3719f39a-021e-42ec-b828-ec22ed6ee6a6.png width="648">