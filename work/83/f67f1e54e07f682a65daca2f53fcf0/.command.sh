#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno170804_leafcutter_merged_sorted.bam -o Pheno170804.junc

grep -v "?" Pheno170804.junc > Pheno170804_Mod.junc

mv Pheno170804_Mod.junc Pheno170804.junc
