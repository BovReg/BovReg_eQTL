#!/bin/bash -ue
stringtie --rf  -p 6 -e -B -G Bos_taurus.ARS-UCD1.2.109.gtf     -o Pheno170907_stringtie.gtf     -A Pheno170907_stringtie.tsv     Pheno170907_merged_sorted.bam
