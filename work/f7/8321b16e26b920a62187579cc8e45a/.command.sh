#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno180108_stringtie.gtf     -A Pheno180108_stringtie.tsv     Pheno180108_merged_sorted.bam
