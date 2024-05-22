#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno171303_stringtie.gtf     -A Pheno171303_stringtie.tsv     Pheno171303_merged_sorted.bam
