# BovReg_eQTL analysis 

## eQTL-Detect-Nextflow based pipeline for eQTL detection 
	- This repository provide the nextflow scripts and demo data to test eQTL analysis and run with large datasets
	- This workflow was developed as a sub-workflow design (in five different nextflow scripts ordered numerically from 00 to 04) to make the workflow distributable across different research units.
	- This workflow was primarily developed for cattle (Bos Taurus), but it can adopted for any other species by changing the reference genome and gtf files
	- User should have a linux environment with the lastest nextflow and docker installed to run these scripts.
	- Care should be taken with the input and output paths while running with new data.
	- The numbering of the scripts is the order in which the scripts should be executed with nextflow.
	- To run a test analysis, demodata is available in the folder with the same name.
	- The script 1,3 and 4 contain json files for calling parameters, whcih can run using the command 
        nextflow run script.nf -params-file file.json 
	- All the software tools required for the pipeline are installed using docker container from docker repositories, some other tools eg: QTLtools can be installed uisng the dockerfile present in the given folder
	- To create the index for the reference genome the fasta and gtf file can be downloaded from
		_fasta: https://ftp.ensembl.org/pub/release-109/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz_
		_gtf: https://ftp.ensembl.org/pub/release-109/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.109.gtf.gz_

## Demo data:
	-The genotype demodata is available in the folder "Demo_genotype_BovReg"
	-The Phenotype data (RNAseq expression data in fastq format) can be downloaded from research data open repository Zenodo  https://zenodo.org/record/7949616
	-The genotype-phenotype corresponding samples can be found in the text file "RNA_WGS_CorresID_BovReg.txt"

## Commands to run scripts:
    -Pre-requisite: Install Docker version >=  20.10.8 and nextflow, Nextflow version > 21.04. 

 	Script 01: nextflow run 01_eQTLDetect.nf -params-file 01_eQTLDetect.json

	Script 02: nextflow run 02_eQTLDetect.nf 

 	Script 03: nextflow run 03_eQTLDetect.nf -params-file 03_eQTLDetect.json

 	Script 04: nextflow run 04_eQTLDetect.nf -params-file 04_eQTLDetect.json

 
