#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno200304_leafcutter_merged_sorted.bam -o Pheno200304.junc

grep -v "?" Pheno200304.junc > Pheno200304_Mod.junc

mv Pheno200304_Mod.junc Pheno200304.junc
