[![Nextflow version](https://img.shields.io/badge/Nextflow-v20.01.0-brightgreen)](https://www.nextflow.io/) [![Docker version](https://img.shields.io/badge/Docker-v20.10.8-blue)](https://www.docker.com/) [![Singularity-ce version](https://img.shields.io/badge/Singularity-v3.11.4-green.svg)](https://www.sylabs.io/) [![Podman version](https://img.shields.io/badge/Podman-v3.4.4-violet.svg)](https://podman.io/)

# eQTL-Detect: Nextflow based pipeline for eQTL detection

eQTL-Detect is a [Nextflow](https://www.nextflow.io/) based bioinformatics workflow to detect the cis, trans and splicing eQTLs (expression quantitative trait loci) by performing associations with genotype and expression data.
This repository provides the required Nextflow-DSL2 scripts to run the pipeline and demo data for trial run. User can run the whole analysis either with a single standalone script or using three separate script modules and all the required tools for running the pipeline can be installed using either docker or singularity or podman container technology. We also provided some of the common configurations for running the pipeline on different high performance clusters (HPC).
This pipeline was primarily developed to detect eQTLs in cattle (Bos taurus), but users can adopt this pipeline for other species by providing the reference genome assembly and transcriptome annotation gtf files of the species of interest.

## Software required
Users need to install [Nextflow](https://www.nextflow.io/) and a container tool, which is either [Docker](https://www.docker.com/) or [Singularity](https://www.sylabs.io/) or [Podman](https://podman.io/).

<img width="500" alt="Fig1" src="https://github.com/user-attachments/assets/4c790fb5-1506-488f-ae9d-a4df8642f0e1"> <img width="600" alt="Fig2" src="https://github.com/user-attachments/assets/e78e92f4-3262-42bc-ad9c-c38594610901">



## Pipeline parameters
The [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config) include all the input parameters to run the pipeline with default values and also the paths for different input files and path for the output directory to store the output results. 
- The input data include the reference genome, reference annotation file and the paths of .tsv files. These .tsv files include the IDs and path of genotype data and expression data..  

##  Input file formats and demo data
Users can download the demo data and can perform a trial run of the pipeline and the links for downloading the test data are given below.
The [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config) include the default parameters along with the paths for the demo data. (NOTE: Also with the input phenotype file user should also declare the type of input files provided with boolean parameters (true or false) based on the available phenotype input (fasta or bam or count matrices) in the [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config). By default the pipeline takes paired-end read FASTA files as input.

- The reference genome and reference annotation file should be provided as fasta format and gtf format respectively. We provided bovine reference genome for trial run which can be downloaded here: [reference genome](https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz) and uncompress using the command "gunzip Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz" and [reference annotation](https://ftp.ensembl.org/pub/release-109/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.109.gtf.gz) and uncompress using the command "gunzip Bos_taurus.ARS-UCD1.2.109.gtf.gz" and save these files in the folder [Reference_Genome](https://github.com/BovReg/BovReg_eQTL/tree/main/Reference_Genome).
  
- The genotype should be provided in [vcf](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format and in a TSV(Tab-Separated Values) text file user should give the chromosome number and file path of the corresponding vcf files.

  For the test run the Demo genotype data can be [downloaded here](https://zenodo.org/records/10997393/files/Demo_genotype_BovReg.tar.gz) and also an examples tsv file [Demodata/Geno_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Geno_input.tsv) was provided. 
  (Note: For the trial run the "Demo_genotype_BovReg.tar.gz" should be saved and uncompressed (tar -xvf Demo_genotype_BovReg.tar.gz) in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/tree/main/Demodata).

- The Phenotype data can be provided in any of the following formats and all the input files paths should be given in a TSV text file. TSV files with demodata of different formats is available in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/tree/main/Demodata).

   1. Raw data (RNAseq expression data in fastq format): Demo_data can be [downloaded here](https://zenodo.org/records/7949616/files/Demo_RNAseqData_BovReg.tar.gz) and a [Demodata/fastq_paired_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/fasta_paired_input.tsv) file was provided with the file paths for different fastq samples present in the demo data. ( Note: For the trial run the downloaded folder containing the corresponding fastq files "Demo_RNAseqData_BovReg.tar.gz" should be saved and uncompressed (tar -xvf Demo_RNAseqData_BovReg.tar.gz)  in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/tree/main/Demodata). 

   2. Aligned reads (RNAseq expression data in bam format): Demo_data can be [downloaded here](https://zenodo.org/records/7950181/files/Demo_RNAseqBam_BovReg.tar.gz) and a [Demodata/Bam_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Bam_input.tsv) file was provided with the file paths for different bam samples present in the demo data. ( Note: For the trial run the downloaded folder containing the corresponding bam files "Demo_RNAseqBam_BovReg.tar.gz" should be saved and uncompressed (tar -xvf Demo_RNAseqBam_BovReg.tar.gz) in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/tree/main/Demodata).

   4. Expression counts across samples for gene level, transcript level and splicing counts (expression count matrices as text file): Demo_data is automatically downloaded during the execution of the [Demodata/Count_matrices.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Count_matrices.tsv) file. 
   
- The genotype-phenotype corresponding samples information should be provided as text file: [Demo_data/RNA_WGS_CorresID_BovReg.txt](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/RNA_WGS_CorresID_BovReg.txt).

## Commands to run the pipeline:

Based on user preferences this analysis can run with a single script or by using modular scripts and they can choose and modify the [config file](https://github.com/BovReg/BovReg_eQTL/tree/main/conf) based on the available computational cluster.

- **Single command approach:** To run the whole pipeline with single command.

        nextflow run main.nf -c conf/env_local.config -profile docker
    
   - The alignment step can be skipped if the user has aligned bam files as input, which can be mentioned as boolean logic 'true' in [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config).

   - This script can also run only by providing the expression count matrices, which can be mentioned as boolean logic 'true' in [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config).

- **Modular approach:** For modular analysis users can opt for the following scripts.
     Note: Users can skip Module 1, if they have aligned and sorted bam files and if users have count matrices the modules 1 and 2 can be ignored and only Module 3 can be used for eQTL detection.

  - Module 1: Indexing the reference genome and aligning the RNAseq reads

        nextflow run module_1_eQTLDetect.nf -c conf/env_local.config -profile docker

  - Module 2: Extract genotypes from samples having corresponding RNAseq data, Quantification and merging RNAseq samples counts.

        nextflow run module_2_eQTLDetect.nf -c conf/env_local.config -profile docker

  - Module 3: Perform cis, trans and sQTL mapping
    
        nextflow run module_3_eQTLDetect.nf -c conf/env_local.config -profile docker

 - **Required parameters**

   - Users should provide the read type and read strandedness for the RNAseq data with boolean logic true or false in the [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config) file.
     - Read type: --pairedEnd_reads, --singleEnd_reads 
     - Strandedness: --firstStranded, --secondStranded and --unStranded
   - Additionally, all the other parameters required to run the pipeline are defined in [nextflow.config](https://github.com/BovReg/BovReg_eQTL/blob/main/nextflow.config) file, users can change the default values if required.
