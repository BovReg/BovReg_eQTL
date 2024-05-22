#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno200304_stringtie.gtf     -A Pheno200304_stringtie.tsv     Pheno200304_merged_sorted.bam
