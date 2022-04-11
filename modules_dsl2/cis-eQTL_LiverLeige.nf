/*
 process cisQTL: The continer used for this process is preperared with the dockerfile
/home/chitneedi/Disk2chitneedi/NextFlow/QTLtools_Dockerfile/Dockerfile
 This dockerfile has the installations for software QTLtools, fatqtl, tabix, vcftools, samtools, bcftools and htslib
 from biocontainers "https://github.com/BioContainers/containers"

 This process performs cis QTL detection across 29 autosome chromsomes using QTL tool
 Nominal pass - to get nominal P-values between your phenotypes and genotypes
 Permutaion test -  in order to get adjusted P-values of association between the phenotypes and the top variants in cis
 The cisQTL detection was perfomed for genes, transcripts and splice variants


*/

process cisQTL_nominal {
 publishDir "${params.outputQTL}/cisQTLResults", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
 tuple val(chr),file (phenotype_g), file (covariate_g), file (genotype_g), file (phenotype_t), file (covariate_t), file (genotype_t) 

 output:

 tuple val(chr), file ('*') , emit:  cis_nominal_qtlResults_ch

 script:

 //output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 //output_s_mod = output_s - ~/FBN_Geno_Chr/
 output_t_mod = output_t - ~/FBN_Geno_Chr/
 output_g_mod = output_g - ~/FBN_Geno_Chr/

 //QTLtools cis --vcf $genotype_e --bed $phenotype_e --cov $covariate_e --nominal 0.01 --out ${output_e}

 // sed '1iPhenotypeID_Exon Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd NominalP-val RegressionSlope Binaryflag' ${output_e} > ${output_e}_cis.nominals.txt
 // rm ${output_e}
 
 //sed -i 's/_liver_Aligned.sortedByCoord.out//g' $phenotype_t
 //sed -i 's/_liver_Aligned.sortedByCoord.out//g' $covariate_t

 """
 QTLtools cis --vcf $genotype_t --bed $phenotype_t --cov $covariate_t --nominal 0.01 --normal --out ${output_t}

 sed '1iPhenotypeID_Transcript Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd NominalP-val RegressionSlope Binaryflag' ${output_t} > FBN_liver_${output_t_mod}_cis.nominal.txt

 rm ${output_t}


 QTLtools cis --vcf $genotype_g --bed $phenotype_g --cov $covariate_g --nominal 0.01 --normal --out ${output_g}

 sed '1iPhenotypeID_Gene Chrmosome start end strand TotalTestVariants Distance VariantId VarinatChrId start end NominalP-val RegressionSlope Binaryflag' ${output_g} >  FBN_liver_${output_g_mod}_cis.nominal.txt

rm ${output_g}

 """
}



/*
 process cisQTL: This process performs cis QTL detection across 29 autosome chromsomes using QTL tool

 The cisQTL detection was perfomed for genes, transcripts and splice variants

*/

process cisQTL_permutation {
 publishDir "${params.outputQTL}/cisQTLResults", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
 tuple val(chr),file (phenotype_g), file (covariate_g), file (genotype_g), file (phenotype_t), file (covariate_t), file (genotype_t) 

 output:

 tuple val(chr), file ('*') , emit:  cis_permu_qtlResults_ch

 script:

 //output_s = genotype_s[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_t = genotype_t[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 output_g = genotype_g[0].toString() - ~/(\.recode)?(\.vcf)?(\.gz)?$/
 //output_s_mod = output_s - ~/FBN_Geno_Chr/
 output_t_mod = output_t - ~/FBN_Geno_Chr/
 output_g_mod = output_g - ~/FBN_Geno_Chr/

 //QTLtools cis --vcf $genotype_e --bed $phenotype_e --cov $covariate_e --permute 10000  --out ${output_e}

 // sed '1iPhenotypeID_Exon Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd NominalP-val RegressionSlope Binaryflag' ${output_e} > ${output_e}_cis.nominals.txt
 // rm ${output_e}
 """

 QTLtools cis --vcf $genotype_t --bed $phenotype_t --cov $covariate_t --normal --permute 1000  --out ${output_t}

 sed '1iPhenotypeID_Transcript Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd DF Dummy 1Parbeta 2Parbeta unknown NominalP-val RegressionSlope EmpiricalP-val adjP-val' ${output_t} > FBN_liver_${output_t_mod}_cis.permu.txt

 rm ${output_t} 


 QTLtools cis --vcf $genotype_g --bed $phenotype_g --cov $covariate_g --normal --permute 1000  --out ${output_g}

 sed '1iPhenotypeID_Gene Chrmosome Start End Strand TotalTestVariants Distance VariantId VarinatChrId VariantStart VariantEnd DF Dummy 1Parbeta 2Parbeta unknown NominalP-val RegressionSlope EmpiricalP-val adjP-val' ${output_g} >  FBN_liver_${output_g_mod}_cis.permu.txt

 rm ${output_g}

 """
}