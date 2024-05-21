[![Nextflow version](https://img.shields.io/badge/Nextflow-v20.01.0-brightgreen)](https://www.nextflow.io/) [![Docker version](https://img.shields.io/badge/Docker-v20.10.8-blue)](https://www.docker.com/) [![Singularity-ce version](https://img.shields.io/badge/Singularity-v3.11.4-green.svg)](https://www.sylabs.io/) [![Podman version](https://img.shields.io/badge/Podman-v3.4.4-violet.svg)](https://podman.io/)


# eQTL-Detect: Nextflow based pipeline for eQTL detection

eQTL-Detect is a Nextflow based bioinformatics workflow to detect the cis, trans and splicing eQTLs (expresssion quantitative trait loci) by perfoming associations with genotype and expression data.
This repository provides the [Nextflow](https://www.nextflow.io/) scripts and demo data to test the eQTL analyses and run with large datasets. In this updated version user can run the whole analysis with a single command or in three separate modules and long with docker the user can opt for either singularity of podman container technologies to install all the required tools for the pipeline.
 It was primarily developed to detect eQTLs in cattle (Bos taurus), but it can be adopted for any other species by providing the reference genome assembly and transcriptome annotation gtf files of the species of interest.

## Software required
- Users need to install  [Nextflow](https://www.nextflow.io/)  and a container tool, which is either [Docker](https://www.docker.com/) or [Singularity](https://www.sylabs.io/) or [Podman](https://podman.io/).


## Required specifications for input files and demo data for trail run
Users can download the demo data and can perform a trail run of the pipeline and the links for downloading the test data are given below.

- The genotype should be provided in [vcf](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format, a demodata for test run can be downloaded here: Demo_genotype_data [download](https://zenodo.org/records/10997393/files/Demo_genotype_BovReg.tar.gz?download=1) and a [Demodata/Geno_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Geno_input.tsv) file was provided with the file paths of different chromosomes.
  
   Note: For the trial run the folder containing all the vcf.gz files "Demo_genotype_BovReg.tar.gz" should be saved and uncompressed in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata). 

- The Phenotype data can be provided in any of the following formats and all the input files paths should be provides in a tsv file. An example tsv files with demodata is available in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/tree/main/Demodata).

   1. Raw data (RNAseq expression data in fastq format): Demo_data [download](https://zenodo.org/records/7949616/files/Demo_RNAseqData_BovReg.tar.gz?download=1) and a [Demodata/fastq_paired_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/fasta_paired_input.tsv) file was provided with the file paths for different fastq samples.
   
      Note: For the trial run the folder containing the corresponding fastq files "Demo_RNAseqData_BovReg.tar.gz" should be saved and uncompressed in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata). 

   2. Aligned reads (RNAseq expression data in bam format): Demo_data [download](https://zenodo.org/records/7950181/files/Demo_RNAseqBam_BovReg.tar.gz?download=1) and a [Demodata/Bam_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Bam_input.tsv) file was provided with the file paths for different bam samples.

      Note: For the trial run the folder containing the corresponding bam files "Demo_RNAseqBam_BovReg.tar.gz" should be saved and uncompressed in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata).

   4. Expression counts across samples for genelevel, transcript level and spilcing counts (expression count matrices as text file): Demo_data [Demodata/Count_matrices.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Count_matrices.tsv). 
   

- The genotype-phenotype corresponding samples information should be provided as text file: [Demo_data/RNA_WGS_CorresID_BovReg.txt](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/RNA_WGS_CorresID_BovReg.txt).

- The reference genome and annotation file for the demo analysis can be downloaded here [reference genome: fasta format ](https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz) and [reference annotation: gtf format](https://ftp.ensembl.org/pub/release-109/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.109.gtf.gz).

## Commands to run scripts:

This analysis can run with a single script or by using modular scripts based on user preferences.

-  Users should provide the read type and read strandedness for the RNAseq data with boolean logic true or false in the [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config) file:

  - Parameters 
    - Read type: --pairedEnd_reads, --singleEnd_reads 
    - Strandedness: --firstStranded, --secondStranded and --unStranded

   **Single command approach:** To run whole pipeline with single command.

   _nextflow run [main.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/main.nf)_ 

- The alignment step can be skipped if the user has aligned bam files as input, which can be mentioned as boolean logic 'true'  in [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config).

- This script can also run only by providing the expression count matrices, which can be mentioned as boolean logic 'true'  in [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config).

  **Modular approach:**

- For modular analysis users can opt for the following scripts. Users can skip Module 1, if they have aligned and sorted bam files and if users have count matrices the modules 1 and 2 can be ignored and only Module 3 can be used for eQTL detection.

  **Module 1:** Indexing the reference genome and aligning the RNAseq reads \
  _nextflow run [module_1_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/module_1_eQTLDetect.nf)_  

  **Module 2:** Extract genotypes from samples having corresponding RNAseq data, Quantification and merging RNAseq samples counts\
  _nextflow run [module_2_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/module_2_eQTLDetect.nf)_

  **Module 3:**  Perform cis, trans and sQTL mapping.\
   _nextflow run [module_3_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/module_3_eQTLDetect.nf)_  

