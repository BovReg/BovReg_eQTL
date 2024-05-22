#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno171310_leafcutter_merged_sorted.bam -o Pheno171310.junc

grep -v "?" Pheno171310.junc > Pheno171310_Mod.junc

mv Pheno171310_Mod.junc Pheno171310.junc
