> Machine Learning and AI Approaches to Multimodal Problems in Computational Biology Hackathon (CMU Libraries and DNAnexus), Mar 03-05, 2025

# Blended_seq_imputation

### Contributors
* Phil Greer
* Tien Ly
* Yuning Zheng
* Gobikrishnan Subramaniam

## Introduction
Blended Genome Exome (BGE) sequencing is an innovative approach developed by the Broad Institute that integrates low-pass whole genome sequencing (WGS) at approximately 3x coverage with 30x coverage whole exome sequencing (WES) on a unified sequencing platform. Unlike traditional genotyping arrays, BGE and low-pass WGS are not limited by predefined probe sets based on specific ancestral data, making them more inclusive for diverse populations. To obtain accurate variant calls for common variants across the genome, both BGE and low-pass WGS data require an imputation step, which relies on high-quality reference panels comprising large, ancestrally diverse samples. In this project, we aim to develop a Nextflow workflow to construct such an imputation reference panel using extensive cohort datasets, leveraging resources like the 1000 Genomes Project (1kGP) and the Human Genome Diversity Project (HGDP), which have been harmonized to create a comprehensive resource of deeply sequenced human genomes.

## Methods
Building upon [GLIMPSE2's tutorial](https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/), which provides bash script snippets for generating reference panels from a single chromosome (specifically, chromosome 22 of the 1000 Genomes Project b38 data from the EBI FTP site), we developed a scalable Nextflow pipeline capable of processing all chromosomes in our dataset. Due to the complexities associated with chromosome X, it was excluded from this project at this time.

Our pipeline utilizes the combined, phased dataset from the HGDP and the 1000 Genomes Project, available in the gnomAD public cloud folders on [Google Cloud 
Platform](gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2) and [AWS](s3://gnomad-public-us-east-1/resources/hgdp_1kg/phased_haplotypes_v2). The dataset is made up of 4091 subjects from diverse ancestries.The pipeline comprises four main steps:

1) Conversion of multiallelic sites: Transform all multiallelic sites into biallelic sites, retaining both single nucleotide polymorphisms (SNPs) and insertions/deletions (indels)
2) Extraction of site information: Extract site information for the entire cohort ignoring specific genotype calls
3) Chunking reference data: Divide the reference data using GLIMPSE2_chunk and prior mapping information provided in the GLIMPSE2 repository
4) Splitting reference chromosomes: Segment the reference chromosomes into binary chunks for all chromosomes

The software tools employed include [bcftools](https://samtools.github.io/bcftools/), [GLIMPSE2](https://odelaneau.github.io/GLIMPSE/), [Nextflow](https://github.com/nextflow-io/nextflow), and [Docker](https://github.com/docker).

### Workflow
![flowchart](figures/nextflow.drawio_whitebackground.png)

## Results
[Slides](https://docs.google.com/presentation/d/1FQdrbxCqVt1jzBpF6MPERIgx8S5RtXtNTYlFa6TiLg4/edit?usp=sharing)

## Operation
Users can execute the Nextflow pipeline on their chosen phased, WGS dataset to create a customized reference panel. The pipeline's modular design facilitates straightforward integration into existing workflows, enhancing its adaptability for various research applications. 

The workflow has been successfully tested on local systems, on-premises HPC systems, DNAnexus, and Google Cloud Platform (GCP). On DNAnexus, the workflow takes < 1 hour to create the final binary reference files for our example dataset using 22 concurrent operators on 4 core virtual machines (mem2_ssd1_v2_x4). The compute cost of generating the reference panel was USD $1.50. 

On DNAnexus, one imports the current git repository similar to the following instrictions for importing [nf-core](https://academy.dnanexus.com/buildingworkflows/nf/importingandrunningnfcore). Once imported, you need a sample sheet with the location of the phased, imputed, genotype data in indexed VCF or BCF format uploaded to your DNAnexus project folder. Please see the samples_dnanexus.txt in the example_inputs folder for proper formatting. You must also upload the map files found in the GLIMPSE2 repository onto your DNAnexus project and provide the workflow with the path to the map files and the output directory for your final reference panel files.  These must also be in the same dx://project-ID:/folder format seen in the samples_dnanexus.txt file. Finally, we reccomend adding "-queue-size 22" to the Nextflow run options window to allow Nextflow to parallelize across all 22 autosomes. You can then run the workflow.

## References

DeFelice M, Grimsby JL, Howrigan D, Yuan K, Chapman SB, Stevens C, DeLuca S, Townsend M, Buxbaum J, Pericak-Vance M, Qin S, Stein DJ, Teferra S, Xavier RJ, Huang H, Martin AR, Neale BM. Blended Genome Exome (BGE) 
as a Cost Efficient Alternative to Deep Whole Genomes or Arrays. bioRxiv [Preprint]. 2024 Jul 30:2024.04.03.587209. doi: 10.1101/2024.04.03.587209. PMID: 38645052; PMCID: PMC11030253.

Koenig Z, Yohannes MT, Nkambule LL, Zhao X, Goodrich JK, Kim HA, Wilson MW, Tiao G, Hao SP, Sahakian N, Chao KR, Walker MA, Lyu Y; gnomAD Project Consortium; Rehm HL, Neale BM, Talkowski ME, Daly MJ, Brand H, 
Karczewski KJ, Atkinson EG, Martin AR. A harmonized public resource of deeply sequenced diverse human genomes. bioRxiv [Preprint]. 2024 Feb 28:2023.01.23.525248. doi: 10.1101/2023.01.23.525248. Update in: Genome 
Res. 2024 Jun 25;34(5):796-809. doi: 10.1101/gr.278378.123. PMID: 36747613; PMCID: PMC9900804.

Rubinacci, S., Ribeiro, D.M., Hofmeister, R.J. et al. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. Nat Genet 53, 120–126 (2021). doi.org/10.1038/s41588-020-00756-0

Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. Nat Genet. 2023 Jul;55(7):1088-1090. doi: 10.1038/s41588-023-01438-3. Epub 2023 
Jun 29. PMID: 37386250; PMCID: PMC10335927.

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 
10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.

Dirk Merkel. 2014. Docker: lightweight Linux containers for consistent development and deployment. Linux J. 2014, 239, Article 2 (March 2014).

