#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno171310_stringtie.gtf     -A Pheno171310_stringtie.tsv     Pheno171310_merged_sorted.bam
