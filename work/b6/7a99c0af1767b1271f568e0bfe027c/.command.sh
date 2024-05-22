#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno211125_leafcutter_merged_sorted.bam -o Pheno211125.junc

grep -v "?" Pheno211125.junc > Pheno211125_Mod.junc

mv Pheno211125_Mod.junc Pheno211125.junc
