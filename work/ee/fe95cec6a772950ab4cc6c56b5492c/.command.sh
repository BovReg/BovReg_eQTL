#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno220410_leafcutter_merged_sorted.bam -o Pheno220410.junc

grep -v "?" Pheno220410.junc > Pheno220410_Mod.junc

mv Pheno220410_Mod.junc Pheno220410.junc
