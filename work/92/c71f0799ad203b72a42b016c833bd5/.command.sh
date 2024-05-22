#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno220116_leafcutter_merged_sorted.bam -o Pheno220116.junc

grep -v "?" Pheno220116.junc > Pheno220116_Mod.junc

mv Pheno220116_Mod.junc Pheno220116.junc
