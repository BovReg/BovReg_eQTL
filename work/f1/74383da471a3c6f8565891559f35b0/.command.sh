#!/bin/bash -ue
## Changes the samples headers in the genotype data similar to RNAseq sampleIDs with correponding sample text file ##

  awk '{print $1}' RNA_WGS_CorresID_BovReg.txt > dataSet_Common_GenoIds.txt

  vcftools --gzvcf Bovreg_demogeno_Chr3.vcf.gz --keep dataSet_Common_GenoIds.txt  --recode  --out dataSet_Geno3

  plink --cow --vcf dataSet_Geno3.recode.vcf  --maf 0.005 --geno 0.5  --mind 0.01 --hwe 1e-6 --keep-allele-order  --recode vcf --out dataSet_QC_Geno3


### rename genotype ids to RNAseq ids ##

  bcftools reheader -s RNA_WGS_CorresID_BovReg.txt -o dataSet_genotypesChr3.vcf    dataSet_QC_Geno3.vcf

  bgzip -f dataSet_genotypesChr3.vcf

  tabix -p vcf  dataSet_genotypesChr3.vcf.gz

   ## Extract the list samples having genotypes to filter count matrices of RNAseq corresponding samples in script 03 #

  bcftools query -l dataSet_genotypesChr3.vcf.gz | awk '{print $1}' > GenoSamples_Chr3.txt
