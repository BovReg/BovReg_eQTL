#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno180108_leafcutter_merged_sorted.bam -o Pheno180108.junc

grep -v "?" Pheno180108.junc > Pheno180108_Mod.junc

mv Pheno180108_Mod.junc Pheno180108.junc
