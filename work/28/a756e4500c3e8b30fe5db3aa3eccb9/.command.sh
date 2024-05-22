#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno220402_stringtie.gtf     -A Pheno220402_stringtie.tsv     Pheno220402_merged_sorted.bam
