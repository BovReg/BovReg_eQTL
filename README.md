# BovReg_eQTL analysis 

## eQTL-Detect: Nextflow based pipeline for eQTL detection 
This repository provide the [Nextflow](https://www.nextflow.io/) scripts and demo data to test eQTL analyses and run with large datasets. It was developed as a sub-workflow design (in five independent nextflow scripts and ordered numerically from 00 to 04) in order make the workflow distributable across different project partners. It was primarily developed for cattle (_Bos taurus_), but it can adopted for any other species by providing the reference genome assembly and transcriptome annotation gtf files of the species of interest. 


## Required software
- [Nextflow](https://www.nextflow.io/) version => 21.04 
- [Docker](https://www.docker.com/) version >=  20.10.8 installed on their machines to run these scripts.

 Care should be taken with the input and output paths of the data. The numbering of the scripts represents the order in which the scripts should be executed. 
                            


To run a test analysis, a demo data is also provided. The script 1,3 and 4 contain json files for calling parameters, which can run using the command _nextflow run script.nf -params-file file.json_. All the software tools required for the pipeline are installed using docker container and for [QTLtools](https://qtltools.github.io/qtltools/) a docker container version was not available, so we generated container using Dockerfile which was also provided in the repository.
The reference genome and annotation file for this demo analysis can be downloaded here [reference genome: fasta format ](https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz) and [reference annotation: gtf format](https://ftp.ensembl.org/pub/release-109/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.109.gtf.gz).
 


  
## Demo data:
- The genotype demodata is available in the folder [Demo_genotype_BovReg](https://github.com/BovReg/BovReg_eQTL/tree/main/Demo_genotype_BovReg). The genotype data can be supplied in raw reads (FastQ), read counts,  and normalized expression levels (e.g., TPM).
- The Phenotype data (RNAseq expression data in fastq format) can be downloaded from research data open repository 
 [Zenodo](https://zenodo.org/record/7949616) 
- The genotype-phenotype corresponding samples can be found in the text file [RNA_WGS_CorresID_BovReg.txt](https://github.com/BovReg/BovReg_eQTL/blob/main/RNA_WGS_CorresID_BovReg.txt)

## Commands to run scripts:

The analysis can run with a single script or using modular scripts based on user preferences.

For single script analysis, the user should use the following command by read type and read strandedness:

- The options include 
   - readtype: --pairedEnd_reads, --singleEnd_reads, 
   - Strandedness: --firstStranded, --secondStranded and --unStranded
**Main Script:** _nextflow run [main.nf] (https://github.com/BovReg/BovReg_eQTL/blob/main/main.nf)
-params-file [main.json] (https://github.com/BovReg/BovReg_eQTL/blob/main/main.json) --pairedEnd_reads --firstStranded_

- If the user has aligned bam files, the alignment step can be skipped using the following command
  - The options include 
     Strandedness: --firstStranded, --secondStranded and --unStranded
**Main Script bam input:** _nextflow run [main.nf] (https://github.com/BovReg/BovReg_eQTL/blob/main/main.nf)
-params-file [main.json] (https://github.com/BovReg/BovReg_eQTL/blob/main/main.json) --bamFiles_input --firstStranded_

- This script can also run by provinding the expression count matrices using the following command
  **Main Script bam input:** _nextflow run [main.nf] (https://github.com/BovReg/BovReg_eQTL/blob/main/main.nf)
-params-file [main.json] (https://github.com/BovReg/BovReg_eQTL/blob/main/main.json) --countMatrices_input_


- For modular analysis users can opt for the follwoing scripts.

  **Script 00:** _nextflow run [00_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/00_eQTLDetect.nf)_ for indexing reference genome.

 **Script 01:** _nextflow run [01_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/01_eQTLDetect.nf) -params-file [01_eQTLDetect.json](https://github.com/BovReg/BovReg_eQTL/blob/main/01_eQTLDetect.json)_ for alignment and quantification of expression data.


**Script 02:**  _nextflow run [02_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/02_eQTLDetect.nf)_ for extracting the genotype data from samples having corresponding RNAseq samples.

 **Script 03:**  _nextflow run [03_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/03_eQTLDetect.nf) -params-file [03_eQTLDetect.json](https://github.com/BovReg/BovReg_eQTL/blob/main/03_eQTLDetect.json)_ for merging the read and transcript counts generated from stringtie and cluster introns found in junction files estimate covriates for splicing sites based on PCs.

 **Script 04:** _nextflow run [04_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/04_eQTLDetect.nf) -params-file [04_eQTLDetect.json](https://github.com/BovReg/BovReg_eQTL/blob/main/04_eQTLDetect.json)_ for performing cis and trans QTL mapping using the RNA seq and corresponding genotype data.
