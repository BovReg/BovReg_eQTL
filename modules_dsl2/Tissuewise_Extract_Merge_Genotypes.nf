

/* 

  If the input data is from the more than one source this process should be used
   - to extract genotype data from different sources

*/


process tissuewise_extractGenotype_Multiple_datasets {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
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



/*

   If the input data is from one source this process should be used
   - to extract genotype data from different sources


*/


process tissuewise_extractGenotype {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(dataSet)
      
      file(sample_dataSet)
      


 output:

      tuple val(chr), file ("dataSet_genotypesChr${chr}.vcf.gz")

      file("*.txt")
      

 script:

 """

 ## Changes the samples headers in the genotype data similar to RNAseq sampleIDs with correponding sample text file ##

  awk '{print \$1}' ${sample_dataSet} > dataSet_Common_GenoIds.txt

  vcftools --gzvcf ${dataSet} --keep dataSet_Common_GenoIds.txt  --recode  --out dataSet_Geno${chr}

  bcftools reheader -s ${sample_dataSet} -o dataSet_genotypesChr${chr}.vcf  \
  dataSet_Geno${chr}.recode.vcf

  bgzip -f dataSet_genotypesChr${chr}.vcf

  tabix -p vcf  dataSet_genotypesChr${chr}.vcf.gz

   ## Extract the list samples having genotypes to filter count matrices of RNAseq corresponding samples in script 03 #

  bcftools query -l dataSet_genotypesChr${chr}.vcf.gz | awk '{print \$1"_stringtie.tsv"}' > GenoSamples_tsv_Chr${chr}.txt

  bcftools query -l dataSet_genotypesChr${chr}.vcf.gz | awk '{print \$1"_stringtie.gtf"}' > GenoSamples_gtf_Chr${chr}.txt

  bcftools query -l dataSet_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamples_splice_Chr${chr}.txt

 """

}




/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=4 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/


process tissuewise_extractGenotype_Rumen {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(FBN), file(LUKE)
      
      file(sample_FBN)

      file(sample_LUKE)
      
   

 output:

      tuple val(chr), file ("Rumen_genotypesChr${chr}.vcf.gz")

      file ("GenoSamplesChr${chr}.txt")
      

 script:

 """


  awk '{print \$1}' ${sample_FBN} > Rumen_FBN_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Rumen_FBN_Common_GenoIds.txt  --recode  --out FBN_Rumen_Geno${chr}

  bcftools reheader -s ${sample_FBN} -o FBN_reheader_Rumen_Geno${chr}.vcf  \
  FBN_Rumen_Geno${chr}.recode.vcf

  bgzip -f FBN_reheader_Rumen_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Rumen_Geno${chr}.vcf.gz


   ### LUKE data set has duplicate RNA seq sample for a genotype sample which has to be collected ##


  cp $sample_LUKE temp

for i in {1..17}
do
  awk '!seen[\$1]++' temp  > Rumen_RNA_WGS_LUKE_CorresID_set\${i}.txt
 

  awk '{print \$1}' Rumen_RNA_WGS_LUKE_CorresID_set\${i}.txt > Rumen_LUKE_Common_GenoIds.txt


  vcftools --vcf ${LUKE} --keep Rumen_LUKE_Common_GenoIds.txt  --recode  --out Rumen_LUKE_Common_GenoIds_set\${i}_Geno${chr}


   bcftools reheader -s Rumen_RNA_WGS_LUKE_CorresID_set\${i}.txt  -o LUKE_reheader_Rumen_Geno${chr}_set\${i}.vcf  \
   Rumen_LUKE_Common_GenoIds_set\${i}_Geno${chr}.recode.vcf

   bgzip -f LUKE_reheader_Rumen_Geno${chr}_set\${i}.vcf 

   tabix -p vcf LUKE_reheader_Rumen_Geno${chr}_set\${i}.vcf.gz


  awk 'NR==FNR {a[\$2]=\$2;next} !(\$2 in a) {print }' Rumen_RNA_WGS_LUKE_CorresID_set\${i}.txt temp > temp2

  mv temp2 temp

 done

  bcftools merge LUKE_reheader_Rumen_Geno${chr}_set*.vcf.gz -o LUKE_reheader_Rumen_Geno${chr}.vcf


  bgzip -f LUKE_reheader_Rumen_Geno${chr}.vcf

  tabix -p vcf  LUKE_reheader_Rumen_Geno${chr}.vcf.gz



   bcftools isec -n=2 -w1  FBN_reheader_Rumen_Geno${chr}.vcf.gz LUKE_reheader_Rumen_Geno${chr}.vcf.gz \
   | bgzip -f -c > Rumen_Merged_Chr${chr}.vcf.gz

  bcftools query -f '%CHROM\t%POS \n' Rumen_Merged_Chr${chr}.vcf.gz > Rumen_CommonVariants${chr}.txt
  
  bcftools filter --regions-file Rumen_CommonVariants${chr}.txt FBN_reheader_Rumen_Geno${chr}.vcf.gz \
  > FBN_Rumen_CommVar_Geno${chr}.vcf

  bgzip -f FBN_Rumen_CommVar_Geno${chr}.vcf

  tabix -p vcf  FBN_Rumen_CommVar_Geno${chr}.vcf.gz

  bcftools filter --regions-file Rumen_CommonVariants${chr}.txt LUKE_reheader_Rumen_Geno${chr}.vcf.gz \
  > LUKE_Rumen_CommVar_Geno${chr}.vcf

  bgzip -f LUKE_Rumen_CommVar_Geno${chr}.vcf

  tabix -p vcf  LUKE_Rumen_CommVar_Geno${chr}.vcf.gz


  bcftools merge -R Rumen_CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > Rumen_genotypesChr${chr}.vcf

  bgzip -f Rumen_genotypesChr${chr}.vcf

  bcftools query -l Rumen_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}



process tissuewise_extractGenotype_Jejunum {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(FBN)
      
      file(sample_FBN)
      
   

 output:

      tuple val(chr), file ("Jejunum_genotypesChr${chr}.vcf.gz")

      file ("GenoSamplesChr${chr}.txt")
      

 script:

 """

  awk '{print \$1}' ${sample_FBN} > Jejunum_FBN_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Jejunum_FBN_Common_GenoIds.txt  --recode  --out FBN_Jejunum_Geno${chr}

  bcftools reheader -s ${sample_FBN} -o Jejunum_genotypesChr${chr}.vcf  \
  FBN_Jejunum_Geno${chr}.recode.vcf

  bgzip -f Jejunum_genotypesChr${chr}.vcf

  tabix -p vcf  Jejunum_genotypesChr${chr}.vcf.gz


  bcftools query -l Jejunum_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}



