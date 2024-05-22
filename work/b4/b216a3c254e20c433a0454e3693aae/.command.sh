#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno220101_leafcutter_merged_sorted.bam -o Pheno220101.junc

grep -v "?" Pheno220101.junc > Pheno220101_Mod.junc

mv Pheno220101_Mod.junc Pheno220101.junc
