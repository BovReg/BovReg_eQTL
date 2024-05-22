#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno211133_leafcutter_merged_sorted.bam -o Pheno211133.junc

grep -v "?" Pheno211133.junc > Pheno211133_Mod.junc

mv Pheno211133_Mod.junc Pheno211133.junc
