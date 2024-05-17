#! /bin/bash

## Bam tsv files ##
ls -l Demodata/Demo_RNAseqBam_BovReg/*_leafcutter_merged_sorted.bam | awk  '{print $9,$9,$9,$9}' | awk '{gsub("Demodata/Demo_RNAseqBam_BovReg/","",$1);gsub("_leafcutter_merged_sorted.bam","",$1);gsub("_leafcutter","",$2);gsub(".bam",".bai",$4); print}' | sed 1i"sampleId stringTieBam leafcutterBam leafcutterBai"  | awk -v OFS="\t" '$1=$1' > Demodata/Bam_input.tsv 


## fastq paired tsv file ##
ls -l Demodata/Demo_RNAseqFastq_BovReg/*R1* | awk '{print $9,$9,$9}' | awk '{gsub("Demodata/Demo_RNAseqFastq_BovReg/","",$1);gsub("_R1_subset.fastq.gz","",$1);gsub("_R1","_R2",$3);print}' | sed 1i"sampleId read1 read2" | awk -v OFS="\t" '$1=$1' > Demodata/fasta_paired_input.tsv

## geno  vcf  files list as tsv file ##
ls -l Demodata/Demo_genotype_BovReg/* | awk '{print $9,$9 }' | awk '{gsub ("Demodata/Demo_genotype_BovReg/Bovreg_demogeno_Chr","",$1);gsub(".vcf.gz","",$1); print }' | sed 1i"chromosome vcfFile" | awk -v OFS="\t" '$1=$1' > Demodata/Geno_input.tsv


