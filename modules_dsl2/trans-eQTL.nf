/*
 process transQTL: This process performs trans QTL detection all the chromsomes using QTL tool
  The transQTL detection was perfomed for genes, transcripts and splice variants

  NOTE: The default htresholds for the parameters --nominal --threshold can be modified based on user requirements.
*/

process trans_eQTL_nominal {
 publishDir "${params.outputQTL}/transQTLResults", mode:'copy'
 container 'praveen/qtltools1.3'

 input:

 tuple val(chr), file (phenotype_g), file (covariate_g),file (genotype_g) ,file (phenotype_t),file (covariate_t),file (genotype_t), file (phenotype_s), file (covariate_s), file (genotype_s)
 
 val(threshold_trans)

 output:

 tuple val(chr), file ('*.best.txt.gz') , emit:  trans_nominal_best_ch

 tuple val(chr), file ('*.bins.txt.gz') , emit:  trans_nominal_bins_ch

 tuple val(chr), file ('*gene_trans.nominal.hits.txt.gz') , emit:  gene_trans_nominal_hits_ch

 tuple val(chr), file ('*transcript_trans.nominal.hits.txt.gz') , emit:  transcript_trans_nominal_hits_ch

 tuple val(chr), file ('*splice_trans.nominal.hits.txt.gz') , emit:  splice_trans_nominal_hits_ch

 script:

 output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_s_mod = output_s - ~/Geno_/
 output_t_mod = output_t - ~/Geno_/
 output_g_mod = output_g - ~/Geno_/

"""

QTLtools trans --vcf $genotype_t --cov $covariate_t --bed $phenotype_t  --normal --nominal  --threshold $threshold_trans --out ${output_t_mod}_trans.nominal

 zcat ${output_t_mod}_trans.nominal.hits.txt.gz | sed '1iTranscript_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_t_mod}_trans.nominal.hits.txt

 gzip -f  ${output_t_mod}_trans.nominal.hits.txt

 zcat ${output_t_mod}_trans.nominal.best.txt.gz |  sed '1iTranscript_Phenotype approxMap P-val'  > ${output_t_mod}_trans.nominal.best.txt

 gzip -f ${output_t_mod}_trans.nominal.best.txt
 

QTLtools trans --vcf $genotype_g --cov $covariate_g --bed $phenotype_g  --normal --nominal  --threshold $threshold_trans --out ${output_g_mod}_trans.nominal

 zcat ${output_g_mod}_trans.nominal.hits.txt.gz | sed '1iGene_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_g_mod}_trans.nominal.hits.txt

 gzip -f  ${output_g_mod}_trans.nominal.hits.txt

 zcat ${output_g_mod}_trans.nominal.best.txt.gz |  sed '1iGene_Phenotype approxMap P-val'  > ${output_g_mod}_trans.nominal.best.txt

 gzip -f ${output_g_mod}_trans.nominal.best.txt

QTLtools trans --vcf $genotype_s --cov $covariate_s --bed $phenotype_s  --normal --nominal  --threshold $threshold_trans --out ${output_s_mod}_trans.nominal

 zcat ${output_s_mod}_trans.nominal.hits.txt.gz | sed '1iTranscript_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_s_mod}_trans.nominal.hits.txt

 gzip -f  ${output_s_mod}_trans.nominal.hits.txt

 zcat ${output_s_mod}_trans.nominal.best.txt.gz |  sed '1iTranscript_Phenotype approxMap P-val'  > ${output_s_mod}_trans.nominal.best.txt

 gzip -f ${output_s_mod}_trans.nominal.best.txt

"""
 
}


/*
 process transQTL: This process performs trans QTL permutations across  all chromsomes using QTL tool
  The transQTL detection was perfomed for genes, transcripts and splice variants

 NOTE: The default htresholds for the parameters --threshold and --permute can be modified based on user requirements.
*/

process trans_eQTL_permu {
 publishDir "${params.outputQTL}/transQTLResults", mode:'copy'
 container 'praveen/qtltools1.3'

 input:


 tuple val(chr),file (phenotype_g), file (covariate_g),file (genotype_g) ,file (phenotype_t),file (covariate_t),file (genotype_t),file (phenotype_s), file (covariate_s), file (genotype_s)
 
 val(permutations_trans)

 val(threshold_trans)
 
 output:

 tuple val(chr), file ('*.best.txt.gz') , emit:  trans_permu_best_ch

 tuple val(chr), file ('*.bins.txt.gz') , emit:  trans_permu_bins_ch

 tuple val(chr), file ('*gene_trans.permu.hits.txt.gz') , emit:  gene_trans_permu_hits_ch

 tuple val(chr), file ('*transcript_trans.permu.hits.txt.gz') , emit:  transcript_trans_permu_hits_ch

 tuple val(chr), file ('*splice_trans.permu.hits.txt.gz') , emit:  splice_trans_permu_hits_ch

 script:

 output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_s_mod = output_s - ~/Geno_/
 output_t_mod = output_t - ~/Geno_/
 output_g_mod = output_g - ~/Geno_/


"""

 QTLtools trans --vcf $genotype_g --bed $phenotype_g --cov $covariate_g  --permute $permutations_trans  --normal --out  ${output_g_mod}_trans.permu


 QTLtools trans --vcf $genotype_t --bed $phenotype_t --cov $covariate_t  --permute $permutations_trans   --normal --out  ${output_t_mod}_trans.permu


 QTLtools trans --vcf $genotype_s --bed $phenotype_s  --cov $covariate_s --permute $permutations_trans  --normal --out  ${output_s_mod}_trans.permu

"""
 
}


/*
 process transQTL: This process performs the FDR estimation of the permuted trans eQTL hits across  all chromsomes using QTL tool
  The transQTL FDR estimation was perfomed for genes, transcripts and splice variants
*/


process trans_eQTL_FDR {
 publishDir "${params.outputQTL}/transQTLResults", mode:'copy'
 container 'praveen/qtltools1.3'

 input:

 tuple val(chr), file(gene_hits), file(gene_hits_permu), file(transcript_hits), file(transcript_hits_permu), file(splice_hits), file(splice_hits_permu)
 
 file(fdr_trans)

 output:

  path ("*.txt") 


 script:


"""

    zcat $gene_hits_permu\
    | sed '1igene_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position Pvalue approxMap Regression_slope' | gzip > gene_hits_permu_mod.gz

    Rscript $fdr_trans $gene_hits gene_hits_permu_mod.gz gene_full_pass_trans_eQTL_FDR_output_Chr${chr}_permu.txt

    sed -i '1igene_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position Pvalue approxMap Regression_slope FDR' gene_full_pass_trans_eQTL_FDR_output_Chr${chr}_permu.txt



    zcat $transcript_hits_permu\
    | sed '1itranscript_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position Pvalue approxMap Regression_slope' | gzip > transcript_hits_permu_mod.gz

    Rscript $fdr_trans $transcript_hits transcript_hits_permu_mod.gz transcript_full_pass_trans_eQTL_FDR_output_Chr${chr}_permu.txt

    sed -i '1itranscript_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position Pvalue approxMap Regression_slope FDR' transcript_full_pass_trans_eQTL_FDR_output_Chr${chr}_permu.txt


    zcat $splice_hits_permu\
    | sed '1isplice_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position Pvalue approxMap Regression_slope' | gzip > splice_hits_permu_mod.gz

    Rscript $fdr_trans $splice_hits splice_hits_permu_mod.gz splice_full_pass_trans_eQTL_FDR_output_Chr${chr}_permu.txt

    sed -i '1isplice_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position Pvalue approxMap Regression_slope FDR' splice_full_pass_trans_eQTL_FDR_output_Chr${chr}_permu.txt


"""
 
}
