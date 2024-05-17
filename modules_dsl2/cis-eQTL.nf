/*
 process cisQTL: The continer used for this process is preperared with the dockerfile
/home/chitneedi/Disk2chitneedi/NextFlow/QTLtools_Dockerfile/Dockerfile
 This dockerfile has the installations for software QTLtools, fatqtl, tabix, vcftools, samtools, bcftools and htslib
 from biocontainers "https://github.com/BioContainers/containers"

 This process performs cis QTL detection across all autosomal chromsomes using QTL tool
 Nominal pass - to get nominal P-values between your phenotypes and genotypes
 Permutaion test -  in order to get adjusted P-values of association between the phenotypes and the top variants in cis
 The cisQTL detection was perfomed for genes, transcripts and splice variants


*/

/* NOTE: The default htresholds for the parameters --nominal can be modified based on user requirements. */

process cisQTL_nominal {
 
 tag "on chromosome ${chr}"
 publishDir "${params.outdir}/cisQTLResults", mode:'copy'
 container 'praveenchitneedi/qtltoolkit:v1.0.1'

 input:
 tuple val(chr), file (phenotype_g), file (covariate_g),file (genotype_g) ,file (phenotype_t),file (covariate_t),file (genotype_t), file (phenotype_s), file (covariate_s), file (genotype_s)
 
 val(nominal_cis)

 output:

 tuple val(chr), file ('*.txt') , emit:  cis_nominal_qtlResults_ch

 script:

 output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_s_mod = output_s - ~/Geno_/
 output_t_mod = output_t - ~/Geno_/
 output_g_mod = output_g - ~/Geno_/



 """
QTLtools cis --bed $phenotype_t --vcf $genotype_t --cov $covariate_t --nominal $nominal_cis --normal --out ${output_t}

sed '1iPhenotypeID_Transcript Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd NominalP-val RegressionSlope Binaryflag' ${output_t} > ${output_t_mod}_cis.nominal.txt

rm ${output_t}


QTLtools cis --bed $phenotype_g --vcf $genotype_g  --cov $covariate_g --nominal $nominal_cis --normal --out ${output_g}

sed '1iPhenotypeID_Gene Chrmosome start end strand TotalTestVariants Distance VariantId VarinatChrId start end NominalP-val RegressionSlope Binaryflag' ${output_g} >  ${output_g_mod}_cis.nominal.txt

rm ${output_g}


 QTLtools cis --bed $phenotype_s --vcf $genotype_s  --cov $covariate_s  --nominal $nominal_cis --normal --out ${output_s}

 sed '1iPhenotypeID_Gene Chrmosome start end strand TotalTestVariants Distance VariantId VarinatChrId start end NominalP-val RegressionSlope Binaryflag' ${output_s} >  ${output_s_mod}_cis.nominal.txt

 rm ${output_s}

"""

 
}



/*
 process cisQTL: This process performs cis QTL detection across all chromsomes using QTL tool

 The cisQTL detection was perfomed for genes, transcripts and splice variants

NOTE: The default number of permutations  --permute 1000 can be modified based on user requirements.

*/

process cisQTL_permutation {
 tag "on chromosome ${chr}"
 publishDir "${params.outdir}/cisQTLResults", mode:'copy'
 container 'praveenchitneedi/qtltoolkit:v1.0.1'

 input:
  tuple val(chr),file (phenotype_g), file (covariate_g),file (genotype_g) ,file (phenotype_t),file (covariate_t),file (genotype_t),file (phenotype_s), file (covariate_s), file (genotype_s)

  val(permutations_cis)

 output:

  tuple val(chr), file ("${output_t_mod}_cis.permu.txt") , emit:  cis_permu_transcritResults_ch

  tuple val(chr), file ("${output_g_mod}_cis.permu.txt") , emit:  cis_permu_geneResults_ch

  tuple val(chr), file ("${output_s_mod}_cis.permu.txt") , emit:  cis_permu_splicingResults_ch

 script:

 output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_s_mod = output_s - ~/Geno_/
 output_t_mod = output_t - ~/Geno_/
 output_g_mod = output_g - ~/Geno_/

"""

 QTLtools cis --vcf $genotype_t  --bed $phenotype_t --cov $covariate_t --normal --permute $permutations_cis --out ${output_t}

 sed '1iPhenotypeID_Transcript Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd DF Dummy 1Parbeta 2Parbeta unknown NominalP-val RegressionSlope EmpiricalP-val adjP-val' ${output_t} > ${output_t_mod}_cis.permu.txt

 rm ${output_t} 


 QTLtools cis --vcf $genotype_g --bed $phenotype_g --cov $covariate_g --normal --permute $permutations_cis  --out ${output_g}

 sed '1iPhenotypeID_Gene Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd DF Dummy 1Parbeta 2Parbeta unknown NominalP-val RegressionSlope EmpiricalP-val adjP-val' ${output_g} >  ${output_g_mod}_cis.permu.txt

 rm ${output_g}


QTLtools cis --vcf $genotype_s  --bed $phenotype_s --cov $covariate_s --normal --permute $permutations_cis  --out ${output_s}

 sed '1iPhenotypeID_Gene Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd DF Dummy 1Parbeta 2Parbeta unknown NominalP-val RegressionSlope EmpiricalP-val adjP-val' ${output_s} >  ${output_s_mod}_cis.permu.txt

 rm ${output_s}
 
"""
 
}






/*
 process cisQTL: This process performs cis QTL detection across all chromosmes using QTL tool

 The cisQTL detection was perfomed for genes, transcripts and splice variants

 NOTE: The default fdr 0.05 can be modified based on user requirements.

*/

process cisQTL_conditional {
 tag "on chromosome ${chr}"
 publishDir "${params.outdir}/cisQTLResults", mode:'copy'
 container 'praveenchitneedi/qtltoolkit:v1.0.1'

 input:
 tuple val(chr),file (phenotype_g), file (covariate_g),file (genotype_g), file(permu_g) ,file (phenotype_t),file (covariate_t),file (genotype_t),file(permu_t),file (phenotype_s), file (covariate_s),file (genotype_s),file(permu_s)

 //file(fdr_cis)

 val(fdr_rate)

 output:

 tuple val(chr), file ("*.txt") 


 script:

 output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_s_mod = output_s - ~/Geno_/
 output_t_mod = output_t - ~/Geno_/
 output_g_mod = output_g - ~/Geno_/


"""


 Rscript $projectDir/bin/runFDR_cis_QTLtools.R  $permu_g 0.05 permu_g_Chr${chr}


 QTLtools cis --vcf $genotype_g --bed $phenotype_g --cov $covariate_g --normal --mapping permu_g_Chr${chr}.thresholds.txt  --out ${output_g}

 sed '1iPheno_gene Chr Start End Strand TotalTestVar Distance VarId VarChrId VarStart VarEnd AssoRank ForwardNomPval ForwardRegSlope Dummy BinaryFlagRank BinFlagPval BackwardNomPval BackRegSlope Dummy BinaryFlagRank BinaryFlagP-val' ${output_g} >  ${output_g_mod}_cis.cond.txt

 rm ${output_g}


 cat ${output_g_mod}_cis.cond.txt| awk '{ if (\$21 == 1) print \$0}' > ${output_g_mod}_cis.cond_topvar.txt 

 sed -i '1iPheno_gene Chr Start End Strand TotalTestVar Distance VarId VarChrId VarStart VarEnd AssoRank ForwardNomPval ForwardRegSlope Dummy BinaryFlagRank BinFlagPval BackwardNomPval BackRegSlope Dummy BinaryFlagRank BinaryFlagP-val' ${output_g_mod}_cis.cond_topvar.txt

 Rscript $projectDir/bin/runFDR_cis_QTLtools.R $permu_t 0.05 permu_t_Chr${chr}


 QTLtools cis --vcf $genotype_t  --bed $phenotype_t --cov $covariate_t  --normal  --mapping permu_t_Chr${chr}.thresholds.txt  --out ${output_t}

 sed '1iPheno_Transcript Chr Start End Strand TotalTestVar Distance VarId VarChrId VarStart VarEnd AssoRank ForwardNomPval ForwardRegSlope Dummy BinaryFlagRank BinFlagPval BackwardNomPval BackRegSlope Dummy BinaryFlagRank BinaryFlagP-val' ${output_t} > ${output_t_mod}_cis.cond.txt

 rm ${output_t} 


 cat ${output_t_mod}_cis.cond.txt| awk '{ if (\$21 == 1) print \$0}' > ${output_t_mod}_cis.cond_topvar.txt

 sed -i '1iPheno_gene Chr Start End Strand TotalTestVar Distance VarId VarChrId VarStart VarEnd AssoRank ForwardNomPval ForwardRegSlope Dummy BinaryFlagRank BinFlagPval BackwardNomPval BackRegSlope Dummy BinaryFlagRank BinaryFlagP-val'  ${output_t_mod}_cis.cond_topvar.txt

 Rscript $projectDir/bin/runFDR_cis_QTLtools.R $permu_s 0.05 permu_s_Chr${chr}

 QTLtools cis --vcf $genotype_s  --bed $phenotype_s --cov $covariate_s  --normal  --mapping permu_s_Chr${chr}.thresholds.txt  --out ${output_s}

 sed '1iPheno_Splice Chr Start End Strand TotalTestVar Distance VarId VarChrId VarStart VarEnd AssoRank ForwardNomPval ForwardRegSlope Dummy BinaryFlagRank BinaFlagPval BackwardNomPval BackRegSlope Dummy BinaryFlagRank BinaryFlagP-val' ${output_s} >  ${output_s_mod}_cis.cond.txt

 rm ${output_s}
 

 cat ${output_s_mod}_cis.cond.txt| awk '{ if (\$21 == 1) print \$0}' > ${output_s_mod}_cis.cond_topvar.txt

 sed -i '1iPheno_gene Chr Start End Strand TotalTestVar Distance VarId VarChrId VarStart VarEnd AssoRank ForwardNomPval ForwardRegSlope Dummy BinaryFlagRank BinFlagPval BackwardNomPval BackRegSlope Dummy BinaryFlagRank BinaryFlagP-val' ${output_s_mod}_cis.cond_topvar.txt

"""
 
}