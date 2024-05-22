#!/bin/bash -ue
regtools junctions extract -s 1 -a 8 -m 50 -M 500000 Pheno170727_leafcutter_merged_sorted.bam -o Pheno170727.junc

grep -v "?" Pheno170727.junc > Pheno170727_Mod.junc

mv Pheno170727_Mod.junc Pheno170727.junc
