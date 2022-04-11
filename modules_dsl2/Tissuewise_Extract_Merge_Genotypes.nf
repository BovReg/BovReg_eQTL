/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=4 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/


process tissuewise_extractGenotype {
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/plink2'

 input:
      tuple val(chr), file(FBN), file(LUKE), file(GIGA), file(UAL)
      
      file(sample_FBN)
      
      file(sample_LUKE)
      
      file(sample_GIGA)
      
      file(sample_UAL)

 output:

      tuple val(chr), file ("Liver_genotypesChr${chr}.vcf.gz"), file("GenoSamplesChr${chr}.txt")
      

 script:

 """

  awk '{print \$1}' ${sample_FBN} > Liver_FBN_Common_GenoIds.txt

  plink2 --vcf ${FBN} --keep Liver_FBN_Common_GenoIds.txt  --cow --recode vcf --out FBN_Liver_Geno${chr}

  bcftools reheader -s ${sample_FBN} -o FBN_reheader_Liver_Geno${chr}.vcf  \
  FBN_Liver_Geno${chr}.vcf

  bgzip -f FBN_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Liver_Geno${chr}.vcf.gz



  awk '{print \$1}' ${sample_LUKE} > Liver_LUKE_Common_GenoIds.txt

  plink2 --vcf ${LUKE} --keep Liver_LUKE_Common_GenoIds.txt  --cow --recode vcf --out LUKE_Liver_Geno${chr}

  bcftools reheader -s ${sample_LUKE}  -o LUKE_reheader_Liver_Geno${chr}.vcf  \
  LUKE_Liver_Geno${chr}.vcf

  bgzip -f LUKE_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  LUKE_reheader_Liver_Geno${chr}.vcf.gz


  awk '{print \$1}' ${sample_GIGA} > Liver_GIGA_Common_GenoIds.txt


  plink2 --vcf ${GIGA} --keep Liver_GIGA_Common_GenoIds.txt --cow --recode vcf --out GIGA_Liver_Geno${chr}

  bcftools reheader -s ${sample_GIGA} -o GIGA_reheader_Liver_Geno${chr}.vcf  GIGA_Liver_Geno${chr}.vcf

  bgzip -f GIGA_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  GIGA_reheader_Liver_Geno${chr}.vcf.gz


  awk '{print \$1}' ${sample_UAL} > Liver_UAL_Changxi_Common_GenoIds.txt


  plink2 --vcf ${UAL} --keep Liver_UAL_Changxi_Common_GenoIds.txt  --cow --recode vcf --out UAL_Changxi_Liver_Geno${chr}

  bcftools reheader -s ${sample_UAL} -o UAL_Changxi_reheader_Liver_Geno${chr}.vcf  UAL_Changxi_Liver_Geno${chr}.vcf

  bgzip -f UAL_Changxi_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  UAL_Changxi_reheader_Liver_Geno${chr}.vcf.gz 

  bcftools isec -n=4 -w1  FBN_reheader_Liver_Geno${chr}.vcf.gz LUKE_reheader_Liver_Geno${chr}.vcf.gz \
  GIGA_reheader_Liver_Geno${chr}.vcf.gz UAL_Changxi_reheader_Liver_Geno${chr}.vcf.gz | bgzip -f -c \
  > Liver_Merged_Chr${chr}.vcf.gz

  bcftools query -f '%CHROM\t%POS \n' Liver_Merged_Chr${chr}.vcf.gz > Liver_CommonVariants${chr}.txt
  
  bcftools filter --regions-file Liver_CommonVariants${chr}.txt FBN_reheader_Liver_Geno${chr}.vcf.gz \
  > FBN_Liver_CommVar_Geno${chr}.vcf

  bgzip -f FBN_Liver_CommVar_Geno${chr}.vcf

  tabix -p vcf  FBN_Liver_CommVar_Geno${chr}.vcf.gz

  bcftools filter --regions-file Liver_CommonVariants${chr}.txt LUKE_reheader_Liver_Geno${chr}.vcf.gz \
  > LUKE_Liver_CommVar_Geno${chr}.vcf

  bgzip -f LUKE_Liver_CommVar_Geno${chr}.vcf

  tabix -p vcf  LUKE_Liver_CommVar_Geno${chr}.vcf.gz

  bcftools filter --regions-file Liver_CommonVariants${chr}.txt GIGA_reheader_Liver_Geno${chr}.vcf.gz \
  > GIGA_Liver_CommVar_Geno${chr}.vcf

  bgzip -f GIGA_Liver_CommVar_Geno${chr}.vcf

  tabix -p vcf  GIGA_Liver_CommVar_Geno${chr}.vcf.gz

  bcftools filter --regions-file Liver_CommonVariants${chr}.txt UAL_Changxi_reheader_Liver_Geno${chr}.vcf.gz \
  > UAL_Changxi_Liver_CommVar_Geno${chr}.vcf

  bgzip -f UAL_Changxi_Liver_CommVar_Geno${chr}.vcf

  tabix -p vcf  UAL_Changxi_Liver_CommVar_Geno${chr}.vcf.gz

  bcftools merge -R Liver_CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > Liver_genotypesChr${chr}.vcf

  bgzip -f Liver_genotypesChr${chr}.vcf

  bcftools query -l Liver_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}