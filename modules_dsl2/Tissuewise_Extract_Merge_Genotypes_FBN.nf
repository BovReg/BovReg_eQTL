/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data

*/

process tissuewise_extractGenotype_Liver_FBN {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(FBN)
      
      file(sample_FBN)
      


 output:

      tuple val(chr), file ("FBN_reheader_Liver_Geno${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """

  awk '{print \$1}' ${sample_FBN} > Liver_FBN_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Liver_FBN_Common_GenoIds.txt  --recode  --out FBN_Liver_Geno${chr}

  bcftools reheader -s ${sample_FBN} -o FBN_reheader_Liver_Geno${chr}.vcf  \
  FBN_Liver_Geno${chr}.recode.vcf

  bgzip -f FBN_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Liver_Geno${chr}.vcf.gz

  bcftools query -l FBN_reheader_Liver_Geno${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt


 """

}