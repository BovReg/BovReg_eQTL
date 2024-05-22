#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno220411_stringtie.gtf     -A Pheno220411_stringtie.tsv     Pheno220411_merged_sorted.bam
