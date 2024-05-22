#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno171012_leafcutter_merged_sorted.bam -o Pheno171012.junc

grep -v "?" Pheno171012.junc > Pheno171012_Mod.junc

mv Pheno171012_Mod.junc Pheno171012.junc
