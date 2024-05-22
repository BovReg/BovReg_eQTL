#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno220116_stringtie.gtf     -A Pheno220116_stringtie.tsv     Pheno220116_merged_sorted.bam
