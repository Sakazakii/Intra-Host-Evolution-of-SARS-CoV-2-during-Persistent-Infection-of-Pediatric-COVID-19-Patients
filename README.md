# **Intra-Host Evolution of SARS-CoV-2 during Persistent Infection of Pediatric COVID-19 Patients**

Charlie R. Boyle, Tien Doan, Estefany Rios-Guzman, Jessica Maciuch, Lacy M. Simons, Dulce Garcia, 
David B. Williams, Arghavan Alisoltani, Egon A. Ozer, Ramon Lorenzo-Redondo, and Judd F. Hultquist
  
## **Table of Contents**
- Introduction
- Significance and Findings
- Study Design
- Dependencies
- Epidemiological Data Sets
- Phylogenetic Data Sets
- Phylogenetic Analysis

## **Introduction**

This repository was curated to store scripts needed to generate both main and supplementary figures in 
Boyle et. al. 2025 (unpublished). The following data are provided:

1. NCBI Accession IDs for sequences used for phylogenetic and phylodynamic analysis
   
2. Codes (.ipynb & .r) used for the analysis of the aforementioned data sets

## **Significance and Findings**
Pediatric patients, especially immunocompromised, exhibit sustained PCR positivity with SARS-CoV-2.

Persistent SARS-CoV-2 positivity generates immune espacape and antiviral resistance non-synonymous mutations.

## **Study Design**
Through a single center retrospective cohort study, we examined the SARS-CoV-2 genetic diversity over a 3-year period in Chicago, Illinois, USA. Between 2020 to 2023, clinical metadata of 20 patient encounters were captured through Ann and Robert H. Lurie Children's Hospital. From March 22nd 2020 to March 22nd 2023, 11,783 residual diagnostic swabs were collected from SARS-CoV-2-positive patients, of which 50 were successfully utilized for SARS-CoV-2 whole-genome sequencing.

## **Dependencies**

### **Python Dependencies**

### **R Dependencies**
- tabplot
- readxl (v.1.4.3)
- dplyr (v.1.1.2)
- lubridate (v.1.9.2)
- ggplot2 v3.5.1
- tidyr (v.1.3.0)
- treeio
- ggtree
- ape
- ggtreeExtra
- ggplot2
- ggsci
- TreeTools (v.1.10.0)
- phylotools (v.1.9-16)
- MAFFT (v7.471)
- scorpio (v0.3.12)
- IQ-Tree (v2.0.5)
- fasttree (v2.1.10)

## **Phylogenetic Analysis**
### **Installing Packages and Activating Environment**
```
module load anaconda3
module load gcc/9.2.0
conda create -c biconda -n dsbr_covid bwa=0.7.17 samtools=1.9 ivar mafft=7.471 fasttree=2.1.10 treetime
conda activate dsbr_covid
mkdir dsbr_covid
cd dsbr_covid
```

### **Concatenate the consensus sequences**
```
cat xxx.fasta > xxx.fa > xxx.fasta
```

### **Alignment**
```
mafft -auto xxx.fasta > xxx.mafft.fasta
```

### **Maximum Likelihood (ML) Phylogenetic Construction (IQTree2)**
```
fasttree -nt -gtr -gamma xxx.mafft.fasta > xxx.tre
```

### **Temporal ML Phylogenetic Construction**
```
treetime --tree xxx.tre --aln xxx.mafft.fasta --dates seqs/dates.csv --outdir treetime
```
