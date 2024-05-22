#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno170804_stringtie.gtf     -A Pheno170804_stringtie.tsv     Pheno170804_merged_sorted.bam
