
/*
   If the input data is from one source this process should be used
   - to extract genotype data from different sources

*/


process tissuewise_extractGenotype {
 tag " on chromosome $chr"
 publishDir "${params.outdir}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(dataSet)
      
      file(sample_dataSet)
      
 output:

      tuple val(chr), file ("dataSet_genotypesChr${chr}.vcf.gz"), emit: genoVcfData_ch

      tuple val(chr), file("GenoSamples_Chr${chr}.txt"), emit: lst_Ind_ch

 script:

 """

 ## Changes the samples headers in the genotype data similar to RNAseq sampleIDs with correponding sample text file ##

  awk '{print \$1}' ${sample_dataSet} > dataSet_Common_GenoIds.txt

  vcftools --gzvcf ${dataSet} --keep dataSet_Common_GenoIds.txt  --recode  --out dataSet_Geno${chr}

### rename genotype ids to RNAseq ids ##

  bcftools reheader -s ${sample_dataSet} -o dataSet_genotypesChr${chr}.vcf  \
  dataSet_Geno${chr}.recode.vcf

  bgzip -f dataSet_genotypesChr${chr}.vcf

  tabix -p vcf  dataSet_genotypesChr${chr}.vcf.gz

   ## Extract the list samples having genotypes to filter count matrices of RNAseq corresponding samples in script 03 #

  bcftools query -l dataSet_genotypesChr${chr}.vcf.gz | awk '{print \$1}' > GenoSamples_Chr${chr}.txt


 """

}






/* 

  If the input data is from the more than one source this process should be used
   - to extract genotype data from different sources

*/

/*
process tissuewise_extractGenotype_Multiple_datasets {
 tag " on chromosome $chr"
 publishDir "${params.outdir}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(dataSet1), file(dataSet2)
      
      file(sample_dataSet1)
      
      file(sample_dataSet2)
      
 output:

      tuple val(chr), file ("dataSet1_dataSet2_genotypesChr${chr}.vcf.gz")

      file("*.txt")
      

 script:

 """
  
  ## Changes the samples headers in the genotype data similar to RNAseq sampleIDs with correponding sample text file ##

  ## Changes the sample IDs in dataset 1 # 

  awk '{print \$1}' ${sample_dataSet1} > dataSet1_Common_GenoIds.txt

  vcftools --gzvcf ${dataSet1} --keep dataSet1_Common_GenoIds.txt  --recode  --out dataSet1_Geno${chr}

  bcftools reheader -s ${sample_dataSet1} -o dataSet1_reheader_Geno${chr}.vcf  \
  dataSet1_Geno${chr}.recode.vcf

  bgzip -f dataSet1_reheader_Geno${chr}.vcf

  tabix -p vcf  dataSet1_reheader_Geno${chr}.vcf.gz


  ## Changes the sample IDs in dataset 2 # 

  awk '{print \$1}' ${sample_dataSet2} > dataSet2_Common_GenoIds.txt

  vcftools --gzvcf ${dataSet2} --keep dataSet2_Common_GenoIds.txt --recode --out dataSet2_Geno${chr}

  bcftools reheader -s ${sample_dataSet2} -o dataSet2_reheader_Geno${chr}.vcf  dataSet2_Geno${chr}.recode.vcf

  bgzip -f dataSet2_reheader_Geno${chr}.vcf

  tabix -p vcf  dataSet2_reheader_Geno${chr}.vcf.gz


  ## NOTE: The above commands can be copied and declare new datasets if the input has more than two different datasets   #

  ## Changes the sample IDs in dataset 3  So on .........#


  ## bcftools isec -n=2 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match # 
  ## NOTE: bcftools isec -n=2 -n=X whrere X signifies number of input datasets # 

  bcftools isec -n=2 -w1  dataSet1_reheader_Geno${chr}.vcf.gz dataSet2_reheader_Geno${chr}.vcf.gz  | bgzip -f -c \
  > Merged_Chr${chr}.vcf.gz

  ## Extracts common variants shared by different datasets # 

  bcftools query -f '%CHROM\t%POS \n' Merged_Chr${chr}.vcf.gz > CommonVariants${chr}.txt

  ## Filters the common vairants from diffrent data sets using bcftools filter comand #
  
  bcftools filter --regions-file CommonVariants${chr}.txt dataSet1_reheader_Geno${chr}.vcf.gz \
  > dataSet1_CommVar_Geno${chr}.vcf

  bgzip -f dataSet1_CommVar_Geno${chr}.vcf

  tabix -p vcf  dataSet1_CommVar_Geno${chr}.vcf.gz


  bcftools filter --regions-file CommonVariants${chr}.txt dataSet2_reheader_Geno${chr}.vcf.gz \
  > dataSet2_CommVar_Geno${chr}.vcf

  bgzip -f dataSet2_CommVar_Geno${chr}.vcf

  tabix -p vcf  dataSet2_CommVar_Geno${chr}.vcf.gz


  ## Merges the datasets after filtering the common variants #

  bcftools merge -R CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > dataSet1_dataSet2_genotypesChr${chr}.vcf

  bgzip -f dataSet1_dataSet2_genotypesChr${chr}.vcf

  
  ## Extract the list samples having genotypes to filter count matrices of RNAseq corresponding samples in script 03 #

  bcftools query -l genotypesChr${chr}.vcf.gz | awk '{print \$1"_stringtie.tsv"}' > GenoSamples_dataSet1_dataSet2_tsv_Chr${chr}.txt

  bcftools query -l genotypesChr${chr}.vcf.gz | awk '{print \$1"_stringtie.gtf"}' > GenoSamples_dataSet1_dataSet2_gtf_Chr${chr}.txt

  bcftools query -l genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamples_dataSet1_dataSet2_splice_Chr${chr}.txt

 """

}

*/




