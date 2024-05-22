#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno220410_stringtie.gtf     -A Pheno220410_stringtie.tsv     Pheno220410_merged_sorted.bam
