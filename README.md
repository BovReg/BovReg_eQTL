# BovReg_eQTL analysis 

## eQTL-Detect: Nextflow based pipeline for eQTL detection 
This repository provide the [Nextflow](https://www.nextflow.io/) scripts and demo data to test eQTL analysis and run with large datasets. It was developed as a sub-workflow design (in five different nextflow scripts ordered numerically from 00 to 04) in order make the workflow distributable across different research units. It was primarily developed for cattle (Bos Taurus), but it can adopted for any other species by giving the reference genome and gtf files of the species of interest. User should have the [Nextflow](https://www.nextflow.io/) version => 21.04 and [Docker](https://www.docker.com/) version >=  20.10.8 installed on their machines to run these scripts. Care should be taken with the input and output paths of the data. The numbering of the scripts represnts the order in which the scripts should be executed. 
                            To run a test analysis, a demo data is available was also provided. The script 1,3 and 4 contain json files for calling parameters, whcih can run using the command _nextflow run script.nf -params-file file.json_. All the software tools required for the pipeline are installed using docker container and for [QTLtools](https://qtltools.github.io/qtltools/) a docker container version was not available, so we generated container using Dockerfile which was also provided in the repository.
The index for the reference genome and annotation file for this demo anlysis can be downloaded here [reference genome: fasta format ](https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz) and [reference annotation: gtf format](https://ftp.ensembl.org/pub/release-109/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.109.gtf.gz).
 
  
## Demo data:
- The genotype demodata is available in the folder [Demo_genotype_BovReg](https://github.com/BovReg/BovReg_eQTL/tree/main/Demo_genotype_BovReg)
- The Phenotype data (RNAseq expression data in fastq format) can be downloaded from research data open repository 
 [Zenodo](https://zenodo.org/record/7949616) 
- The genotype-phenotype corresponding samples can be found in the text file [RNA_WGS_CorresID_BovReg.txt](https://github.com/BovReg/BovReg_eQTL/blob/main/RNA_WGS_CorresID_BovReg.txt)

## Commands to run scripts:
 **Pre-requisite:** Install [Docker](https://www.docker.com/) version >=  20.10.8 and  [Nextflow](https://www.nextflow.io/) version => 21.04. 

  **Script 00:** _nextflow run [00_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/00_eQTLDetect.nf)_

 **Script 01:** _nextflow run [01_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/01_eQTLDetect.nf) -params-file [01_eQTLDetect.json](https://github.com/BovReg/BovReg_eQTL/blob/main/01_eQTLDetect.json)_

**Script 02:**  _nextflow run [02_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/02_eQTLDetect.nf)_ 

 **Script 03:**  _nextflow run [03_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/03_eQTLDetect.nf) -params-file [03_eQTLDetect.json](https://github.com/BovReg/BovReg_eQTL/blob/main/03_eQTLDetect.json)_

 **Script 04:** _nextflow run [04_eQTLDetect.nf](https://github.com/BovReg/BovReg_eQTL/blob/main/04_eQTLDetect.nf) -params-file [04_eQTLDetect.json](https://github.com/BovReg/BovReg_eQTL/blob/main/04_eQTLDetect.json)_

 
