#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno210103_stringtie.gtf     -A Pheno210103_stringtie.tsv     Pheno210103_merged_sorted.bam
