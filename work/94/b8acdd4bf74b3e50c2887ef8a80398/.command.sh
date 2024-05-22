#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno220411_leafcutter_merged_sorted.bam -o Pheno220411.junc

grep -v "?" Pheno220411.junc > Pheno220411_Mod.junc

mv Pheno220411_Mod.junc Pheno220411.junc
