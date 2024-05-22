#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno170727_stringtie.gtf     -A Pheno170727_stringtie.tsv     Pheno170727_merged_sorted.bam
