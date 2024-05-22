#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno220112_stringtie.gtf     -A Pheno220112_stringtie.tsv     Pheno220112_merged_sorted.bam
