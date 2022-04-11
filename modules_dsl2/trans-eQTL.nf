/*
 process transQTL: This process performs trans QTL detection across 29 autosome chromsomes using QTL tool
  The transQTL detection was perfomed for genes, transcripts and splice variants
*/

process transQTL {
 publishDir "${params.outputQTL}/transQTLResults", mode:'copy'
 container 'praveen/qtltools'

 input:

 file phenotype_g from trans_phenotype_gene_Bed_ch

 file genotype_g from trans_genotypeQTL_gene_ch

 file phenotype_t from trans_phenotype_transcript_Bed_ch

 file genotype_t from trans_genotypeQTL_transcript_ch

 //file phenotype_t from trans_phenotype_exon_Bed_ch

 //file genotype_t from trans_genotype_exon_QTL_ch

 file phenotype_s from trans_phenotype_splice_Bed_ch

 file genotype_s from trans_genotypeQTL_splice_ch

 output:

 file '*' , emit:  trans_qtlResults_ch

 script:

 output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_s_mod = output_s - ~/FBN_Geno_/
 output_t_mod = output_t - ~/FBN_Geno_/
 output_g_mod = output_g - ~/FBN_Geno_/

 //QTLtools trans --vcf $genotype_e --bed $phenotype_e  --nominal 0.01 --threshold 1e-5 --out ${output_e}_trans_exon.nominals

 //zcat ${output_e}_trans_exon.nominals.hits.txt.gz | sed '1iExon_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_e}_trans_exon.nominals.hits.txt

 //gzip -f ${output_e}_trans_exon.nominals.hits.txt

 //zcat ${output_e}_trans_exon.nominals.best.txt.gz |  sed '1iExon_Phenotype approxMap P-val'  > ${output_e}_trans_exon.nominals.best.txt

 //gzip -f ${output_e}_trans_exon.nominals.best.txt
 """
 QTLtools trans --vcf $genotype_t --bed $phenotype_t  --nominal 0.01 --threshold 1e-5 --out ${output_t_mod}_trans.nominals

 zcat ${output_t_mod}_trans.nominals.hits.txt.gz | sed '1iTranscript_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_t_mod}_trans.nominals.hits.txt

 gzip -f  ${output_t_mod}_trans.nominals.hits.txt

 zcat ${output_t_mod}_trans.nominals.best.txt.gz |  sed '1iTranscript_Phenotype approxMap P-val'  > ${output_t_mod}_trans.nominals.best.txt

 gzip -f ${output_t_mod}_trans.nominals.best.txt
 

 QTLtools trans --vcf $genotype_g --bed $phenotype_g  --nominal 0.01 --threshold 1e-5 --out ${output_g_mod}_trans.nominals

 zcat ${output_g_mod}_trans.nominals.hits.txt.gz | sed '1iGene_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_g_mod}_trans.nominals.hits.txt

 gzip -f  ${output_g_mod}_trans.nominals.hits.txt

 zcat ${output_g_mod}_trans.nominals.best.txt.gz |  sed '1iGene_Phenotype approxMap P-val'  > ${output_g_mod}_trans.nominals.best.txt

 gzip -f ${output_g_mod}_trans.nominals.best.txt


 QTLtools trans --vcf $genotype_t --bed $phenotype_t  --nominal 0.01 --threshold 1e-5 --out ${output_s_mod}_trans.nominals

 zcat ${output_s_mod}_trans.nominals.hits.txt.gz | sed '1iTranscript_Phenotype Phenotype_Chr Phenotype_Start Varinat_ID Variant_Chr Variant_Position P-value approxMap Regression_slope' > ${output_s_mod}_trans.nominals.hits.txt

 gzip -f  ${output_s_mod}_trans.nominals.hits.txt

 zcat ${output_s_mod}_trans.nominals.best.txt.gz |  sed '1iTranscript_Phenotype approxMap P-val'  > ${output_s_mod}_trans.nominals.best.txt

 gzip -f ${output_s_mod}_trans.nominals.best.txt

 """
}
