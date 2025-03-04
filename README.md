# Blended_seq_imputation
### Contributors
* Phil Greer
* Kara 
* Tien Ly
* Yuning Zheng

## Introduction

Blended Genome Exome sequencing is a new method developed by the Broad to aquire both low pass whole genome sequencing data (~3x depth) with 30x 
Whole exome sequencing data on a combined sequencing platform. Unlike array data, BGE and low pass WGS are ancestry agnostic in that the 
variants information from the sample is not a pre determined set of probes generated fro prior ancestry specific data. Both low pass WGS and BGE 
sequencing data requires an imputation step to get reliable variant calls for commoon variants across the whole genome. Imputation is dependent on 
high quality reference panels consisting of large samples of diverse ancestry. In this project, we will generate a nextflow workflow to generate an 
imputation reference panel on large scale cohort files. 

## Methods
Phased WGS data that includes multiple ancestries is needed to build an adequate reference panel. 
Using the Human Genome Diversity Project + 1000 Genomes combined, phased dataset, we have implemented a nextflow pipeline consisting of 4 steps. 

1) Convert all multiallelic site to biallelic sites, keeping both SNPs and indels.

2) extract the site information for the entire cohort. 

3) chunk out the reference data using GLIPMSE2_chunk. 

4) split out the reference chromosomes into binary chuks for all chromosomes.

Possible testing based on available BGE data

Software used includes bcftools (https://samtools.github.io/bcftools/) and Glimpse2 (https://odelaneau.github.io/GLIMPSE/)

### Data
Four our inital test. we will be using the Human Genome Diversity Project + 1000 Genomes combined dataset. The phased version of the dataset can be 
found on the gnomAD publis cloud folders on Google: gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2 and on AWS: 
s3://gnomad-public-us-east-1/resources/hgdp_1kg/phased_haplotypes_v2

### Workflow
![flowchart](figures/flowchart.png)
## Results

## Discussion

## References

