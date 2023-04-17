/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=4 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/

/*
process tissuewise_extractGenotype_Liver {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(FBN), file(LUKE), file(GIGA), file(UAL)
      
      file(sample_FBN)
      
      file(sample_LUKE)
      
      file(sample_GIGA)
      
      file(sample_UAL)

 output:

      tuple val(chr), file ("Liver_genotypesChr${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """

  awk '{print \$1}' ${sample_FBN} > Liver_FBN_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Liver_FBN_Common_GenoIds.txt  --recode  --out FBN_Liver_Geno${chr}

  bcftools reheader -s ${sample_FBN} -o FBN_reheader_Liver_Geno${chr}.vcf  \
  FBN_Liver_Geno${chr}.recode.vcf

  bgzip -f FBN_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Liver_Geno${chr}.vcf.gz




  ### LUKE data set has duplicate RNA seq sample for a genotype sample which has to be collected ##


  cp $sample_LUKE temp

for i in {1..17}
do
  awk '!seen[\$1]++' temp  > Liver_RNA_WGS_LUKE_CorresID_set\${i}.txt
 

  awk '{print \$1}' Liver_RNA_WGS_LUKE_CorresID_set\${i}.txt > Liver_LUKE_Common_GenoIds.txt


  vcftools --vcf ${LUKE} --keep Liver_LUKE_Common_GenoIds.txt  --recode  --out Liver_LUKE_Common_GenoIds_set\${i}_Geno${chr}


   bcftools reheader -s Liver_RNA_WGS_LUKE_CorresID_set\${i}.txt  -o LUKE_reheader_Liver_Geno${chr}_set\${i}.vcf  \
   Liver_LUKE_Common_GenoIds_set\${i}_Geno${chr}.recode.vcf

   bgzip -f LUKE_reheader_Liver_Geno${chr}_set\${i}.vcf 

   tabix -p vcf LUKE_reheader_Liver_Geno${chr}_set\${i}.vcf.gz


  awk 'NR==FNR {a[\$2]=\$2;next} !(\$2 in a) {print }' Liver_RNA_WGS_LUKE_CorresID_set\${i}.txt temp > temp2

  mv temp2 temp

 done



  bcftools merge LUKE_reheader_Liver_Geno${chr}_set*.vcf.gz -o LUKE_reheader_Liver_Geno${chr}.vcf


  bgzip -f LUKE_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  LUKE_reheader_Liver_Geno${chr}.vcf.gz





  awk '{print \$1}' ${sample_GIGA} > Liver_GIGA_Common_GenoIds.txt


  vcftools --gzvcf ${GIGA} --keep Liver_GIGA_Common_GenoIds.txt --recode --out GIGA_Liver_Geno${chr}

  bcftools reheader -s ${sample_GIGA} -o GIGA_reheader_Liver_Geno${chr}.vcf  GIGA_Liver_Geno${chr}.recode.vcf

  bgzip -f GIGA_reheader_Liver_Geno${chr}.vcf

  tabix -p vcf  GIGA_reheader_Liver_Geno${chr}.vcf.gz


  awk '{print \$1}' ${sample_UAL} > Liver_UAL_Changxi_Common_GenoIds.txt


  vcftools --gzvcf ${UAL} --keep Liver_UAL_Changxi_Common_GenoIds.txt  --recode --out UAL_Changxi_Liver_Geno${chr}

  bcftools reheader -s ${sample_UAL} -o UAL_Changxi_reheader_Liver_Geno${chr}.vcf  UAL_Changxi_Liver_Geno${chr}.recode.vcf

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

*/

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



/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=3 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/


process tissuewise_extractGenotype_Muscle {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
      tuple val(chr), file(FBN), file(UAL)
      
      file(sample_FBN)
      
      file(sample_UAL_LL)

      file(sample_UAL_GM)

 output:

      tuple val(chr), file ("Muscle_genotypesChr${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """

  awk '{print \$1}' ${sample_FBN} > Muscle_FBN_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Muscle_FBN_Common_GenoIds.txt  --recode  --out FBN_Muscle_Geno${chr}

  bcftools reheader -s ${sample_FBN} -o FBN_reheader_Muscle_Geno${chr}.vcf  \
  FBN_Muscle_Geno${chr}.recode.vcf

  bgzip -f FBN_reheader_Muscle_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Muscle_Geno${chr}.vcf.gz

  ## +++NOTE: two UAL samples K021B and K047B were named differently in genotype data (K201400021, K201400047) were considered 

  awk '{print \$1}' ${sample_UAL_LL} > Muscle_UAL_Changxi_Common_GenoIds.txt

  vcftools --gzvcf ${UAL} --keep Muscle_UAL_Changxi_Common_GenoIds.txt   --recode  --out UAL_Changxi_Muscle_Geno${chr}

  bcftools reheader -s ${sample_UAL_LL} -o UAL_Changxi_LL_reheader_Muscle_Geno${chr}.vcf  UAL_Changxi_Muscle_Geno${chr}.recode.vcf

  bgzip -f UAL_Changxi_LL_reheader_Muscle_Geno${chr}.vcf

  tabix -p vcf  UAL_Changxi_LL_reheader_Muscle_Geno${chr}.vcf.gz 


  awk '{print \$1}' ${sample_UAL_GM} > Muscle_UAL_Changxi_Common_GenoIds.txt

  vcftools --gzvcf ${UAL} --keep Muscle_UAL_Changxi_Common_GenoIds.txt  --recode  --out UAL_Changxi_Muscle_Geno${chr}

  bcftools reheader -s ${sample_UAL_GM} -o UAL_Changxi_GM_reheader_Muscle_Geno${chr}.vcf  UAL_Changxi_Muscle_Geno${chr}.recode.vcf

  bgzip -f UAL_Changxi_GM_reheader_Muscle_Geno${chr}.vcf

  tabix -p vcf  UAL_Changxi_GM_reheader_Muscle_Geno${chr}.vcf.gz 

  ### --- Take care of isec -n when having the intersection of different vcf files ---

  bcftools isec -n=3 -w1  FBN_reheader_Muscle_Geno${chr}.vcf.gz UAL_Changxi_LL_reheader_Muscle_Geno${chr}.vcf.gz \
  UAL_Changxi_GM_reheader_Muscle_Geno${chr}.vcf.gz  | bgzip -f -c \
  > Muscle_Merged_Chr${chr}.vcf.gz

  bcftools query -f '%CHROM\t%POS \n' Muscle_Merged_Chr${chr}.vcf.gz > Muscle_CommonVariants${chr}.txt
  
  bcftools filter --regions-file Muscle_CommonVariants${chr}.txt FBN_reheader_Muscle_Geno${chr}.vcf.gz \
  > FBN_Muscle_CommVar_Geno${chr}.vcf

  bgzip -f FBN_Muscle_CommVar_Geno${chr}.vcf

  tabix -p vcf  FBN_Muscle_CommVar_Geno${chr}.vcf.gz


  bcftools filter --regions-file Muscle_CommonVariants${chr}.txt UAL_Changxi_LL_reheader_Muscle_Geno${chr}.vcf.gz \
  > UAL_Changxi_LL_Muscle_CommVar_Geno${chr}.vcf

  bgzip -f UAL_Changxi_LL_Muscle_CommVar_Geno${chr}.vcf

  tabix -p vcf UAL_Changxi_LL_Muscle_CommVar_Geno${chr}.vcf.gz


  bcftools filter --regions-file Muscle_CommonVariants${chr}.txt UAL_Changxi_GM_reheader_Muscle_Geno${chr}.vcf.gz \
  > UAL_Changxi_GM_Muscle_CommVar_Geno${chr}.vcf

  bgzip -f UAL_Changxi_GM_Muscle_CommVar_Geno${chr}.vcf

  tabix -p vcf  UAL_Changxi_GM_Muscle_CommVar_Geno${chr}.vcf.gz


  bcftools merge -R Muscle_CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > Muscle_genotypesChr${chr}.vcf

  bgzip -f Muscle_genotypesChr${chr}.vcf

  bcftools query -l Muscle_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}




/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=2 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/


process tissuewise_extractGenotype_Adipose {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools'

 input:
      tuple val(chr), file(LUKE), file(UAL)
      
      file(sample_LUKE)
      
      file(sample_UAL)

   

 output:

      tuple val(chr), file ("Adipose_genotypesChr${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """

  #awk '{print \$1}' ${sample_LUKE} > Adipose_LUKE_Common_GenoIds.txt

  ### remove duplicate variants using bcftools norm -d both ###

  #bcftools norm -d both ${LUKE} -Oz -o dedup_${LUKE}

  #vcftools --gzvcf dedup_${LUKE} --keep Adipose_LUKE_Common_GenoIds.txt  --recode  --out LUKE_Adipose_Geno${chr}

  #bcftools reheader -s ${sample_LUKE} -o LUKE_reheader_Adipose_Geno${chr}.vcf  \
  #LUKE_Adipose_Geno${chr}.recode.vcf

  #bgzip -f LUKE_reheader_Adipose_Geno${chr}.vcf

  #tabix -p vcf  LUKE_reheader_Adipose_Geno${chr}.vcf.gz




  ### LUKE data set has duplicate RNA seq sample for a genotype sample which has to be collected ##


  cp $sample_LUKE temp

for i in {1..17}
do
  awk '!seen[\$1]++' temp  > Adipose_RNA_WGS_LUKE_CorresID_set\${i}.txt
 

  awk '{print \$1}' Adipose_RNA_WGS_LUKE_CorresID_set\${i}.txt > Adipose_LUKE_Common_GenoIds.txt


  vcftools --vcf ${LUKE} --keep Adipose_LUKE_Common_GenoIds.txt  --recode  --out Adipose_LUKE_Common_GenoIds_set\${i}_Geno${chr}


   bcftools reheader -s Adipose_RNA_WGS_LUKE_CorresID_set\${i}.txt  -o LUKE_reheader_Adipose_Geno${chr}_set\${i}.vcf  \
   Adipose_LUKE_Common_GenoIds_set\${i}_Geno${chr}.recode.vcf

   bgzip -f LUKE_reheader_Adipose_Geno${chr}_set\${i}.vcf 

   tabix -p vcf LUKE_reheader_Adipose_Geno${chr}_set\${i}.vcf.gz


  awk 'NR==FNR {a[\$2]=\$2;next} !(\$2 in a) {print }' Adipose_RNA_WGS_LUKE_CorresID_set\${i}.txt temp > temp2

  mv temp2 temp

 done



  bcftools merge LUKE_reheader_Adipose_Geno${chr}_set*.vcf.gz -o LUKE_reheader_Adipose_Geno${chr}.vcf


  bgzip -f LUKE_reheader_Adipose_Geno${chr}.vcf

  tabix -p vcf  LUKE_reheader_Adipose_Geno${chr}.vcf.gz








  ## +++NOTE: two UAL samples K021B and K047B were named differently in genotype data (K201400021, K201400047) were considered 

  awk '{print \$1}' ${sample_UAL} > Adipose_UAL_Changxi_Common_GenoIds.txt

  ### remove duplicate variants using bcftools norm -d both ###
  bcftools norm -d both ${UAL} -Oz -o dedup_${UAL}

  vcftools --gzvcf dedup_${UAL} --keep Adipose_UAL_Changxi_Common_GenoIds.txt   --recode  --out UAL_Changxi_Adipose_Geno${chr}

  bcftools reheader -s ${sample_UAL} -o UAL_Changxi_reheader_Adipose_Geno${chr}.vcf  UAL_Changxi_Adipose_Geno${chr}.recode.vcf

  bgzip -f UAL_Changxi_reheader_Adipose_Geno${chr}.vcf

  tabix -p vcf  UAL_Changxi_reheader_Adipose_Geno${chr}.vcf.gz 



  

  ### --- Take care of isec -n when having the intersection of different vcf files ---

  bcftools isec -n=2 -w1  LUKE_reheader_Adipose_Geno${chr}.vcf.gz UAL_Changxi_reheader_Adipose_Geno${chr}.vcf.gz \
  | bgzip -f -c > Adipose_Merged_Chr${chr}.vcf.gz

  bcftools query -f '%CHROM\t%POS \n' Adipose_Merged_Chr${chr}.vcf.gz > Adipose_CommonVariants${chr}.txt
  
  bcftools filter  --targets-file Adipose_CommonVariants${chr}.txt LUKE_reheader_Adipose_Geno${chr}.vcf.gz \
  > LUKE_Adipose_CommVar_Geno${chr}.vcf

  ## bcftools norm -d both LUKE_Adipose_CommVar_Geno${chr}.vcf  -o temp

  ## mv temp LUKE_Adipose_CommVar_Geno${chr}.vcf

  bgzip -f LUKE_Adipose_CommVar_Geno${chr}.vcf

  tabix -p vcf  LUKE_Adipose_CommVar_Geno${chr}.vcf.gz

 
  ### -R option takes into account overlapping records. If a strict subset by position is required, add (or replace with) the -T option

  bcftools filter  --targets-file Adipose_CommonVariants${chr}.txt UAL_Changxi_reheader_Adipose_Geno${chr}.vcf.gz \
  > UAL_Changxi_Adipose_CommVar_Geno${chr}.vcf

  ## bcftools norm -d both UAL_Changxi_Adipose_CommVar_Geno${chr}.vcf  -o temp

  ## mv temp UAL_Changxi_Adipose_CommVar_Geno${chr}.vcf

  bgzip -f UAL_Changxi_Adipose_CommVar_Geno${chr}.vcf

  tabix -p vcf UAL_Changxi_Adipose_CommVar_Geno${chr}.vcf.gz



  ## bcftools merge -R Adipose_CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > Adipose_genotypesChr${chr}.vcf

  bcftools merge *CommVar_Geno${chr}.vcf.gz > Adipose_genotypesChr${chr}.vcf


  bgzip -f Adipose_genotypesChr${chr}.vcf

  bcftools query -l Adipose_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}



/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=2 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/


process tissuewise_extractGenotype_Blood {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools'

 input:
      tuple val(chr), file(koch_FBN),file(FBN),  file(GIGA), file(UAL), file(ag_vic)
      
      file(sample_FBN)
      
      file(sample_FBNkoch)

      file(sample_GIGA)

      file(sample_UAL)

      file(sample_agvic)

   

 output:

      tuple val(chr), file ("Blood_genotypesChr${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """
  ### Koch data set has duplicate RNA seq sample for a genotype sample which has to be collected ##

  awk '{print \$1}' ${sample_FBNkoch} > Blood_koch_FBN_Common_GenoIds.txt

  awk '!a[\$1]++' ${sample_FBNkoch} > set1_kocksample.txt

  awk 'NR==FNR{L[\$1]=FNR; next} L[\$1]==FNR'  ${sample_FBNkoch} ${sample_FBNkoch} \
   > set2_kocksample.txt

  vcftools --gzvcf ${koch_FBN} --keep Blood_koch_FBN_Common_GenoIds.txt  --recode  --out koch_FBN_Blood_Geno${chr}

  bcftools reheader -s set1_kocksample.txt  -o koch_FBN_reheader_Blood_Geno${chr}_set1.vcf  \
  koch_FBN_Blood_Geno${chr}.recode.vcf

  bgzip -f koch_FBN_reheader_Blood_Geno${chr}_set1.vcf 

  tabix -p vcf koch_FBN_reheader_Blood_Geno${chr}_set1.vcf.gz



  bcftools reheader -s set2_kocksample.txt  -o koch_FBN_reheader_Blood_Geno${chr}_set2.vcf  \
  koch_FBN_Blood_Geno${chr}.recode.vcf


   bgzip -f koch_FBN_reheader_Blood_Geno${chr}_set2.vcf 

  tabix -p vcf koch_FBN_reheader_Blood_Geno${chr}_set2.vcf.gz


  bcftools merge koch_FBN_reheader_Blood_Geno${chr}_set1.vcf.gz koch_FBN_reheader_Blood_Geno${chr}_set2.vcf.gz \
  -o koch_FBN_reheader_Blood_Geno${chr}.vcf


  bgzip -f koch_FBN_reheader_Blood_Geno${chr}.vcf

  tabix -p vcf  koch_FBN_reheader_Blood_Geno${chr}.vcf.gz



  ### FBN data set has three time with same genotype 00h, 24h, 96h duplicate RNA seq sample for a genotype sample which has to be collected ##

  grep '00h' ${sample_FBN} > Blood_RNA_WGS_FBN_CorresID_00h.txt
  
  grep '24h' ${sample_FBN} > Blood_RNA_WGS_FBN_CorresID_24h.txt

  grep '96h' ${sample_FBN} > Blood_RNA_WGS_FBN_CorresID_96h.txt

  awk '{print \$1}' ${sample_FBN} > Blood_FBN_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Blood_FBN_Common_GenoIds.txt  --recode  --out FBN_Blood_Geno${chr}

 
 awk '{print \$2}' Blood_RNA_WGS_FBN_CorresID_00h.txt > Blood_RNA_IDs_00h.txt

 awk '{print \$2}' Blood_RNA_WGS_FBN_CorresID_24h.txt > Blood_RNA_IDs_24h.txt

 awk '{print \$2}' Blood_RNA_WGS_FBN_CorresID_96h.txt > Blood_RNA_IDs_96h.txt

  bcftools reheader -s Blood_RNA_WGS_FBN_CorresID_00h.txt -o FBN_reheader_Blood_00h_Geno${chr}.vcf  \
  FBN_Blood_Geno${chr}.recode.vcf

  vcftools --vcf FBN_reheader_Blood_00h_Geno${chr}.vcf --keep Blood_RNA_IDs_00h.txt --recode \
   --out temp_00h_${chr}

  bgzip -f  temp_00h_${chr}.recode.vcf

  tabix -p vcf temp_00h_${chr}.recode.vcf.gz

  bcftools reheader -s Blood_RNA_WGS_FBN_CorresID_24h.txt -o FBN_reheader_Blood_24h_Geno${chr}.vcf  \
  FBN_Blood_Geno${chr}.recode.vcf

  vcftools --vcf FBN_reheader_Blood_24h_Geno${chr}.vcf --keep Blood_RNA_IDs_24h.txt --recode \
   --out temp_24h_${chr}

  bgzip -f  temp_24h_${chr}.recode.vcf

  tabix -p vcf temp_24h_${chr}.recode.vcf.gz

 bcftools reheader -s Blood_RNA_WGS_FBN_CorresID_96h.txt -o FBN_reheader_Blood_96h_Geno${chr}.vcf  \
  FBN_Blood_Geno${chr}.recode.vcf

  vcftools --vcf FBN_reheader_Blood_96h_Geno${chr}.vcf --keep Blood_RNA_IDs_96h.txt --recode \
   --out temp_96h_${chr}

  bgzip -f  temp_96h_${chr}.recode.vcf

  tabix -p vcf temp_96h_${chr}.recode.vcf.gz

  
  bcftools merge temp_00h_${chr}.recode.vcf.gz temp_24h_${chr}.recode.vcf.gz temp_96h_${chr}.recode.vcf.gz \
  -o FBN_reheader_Blood_Geno${chr}.vcf

  bgzip -f FBN_reheader_Blood_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Blood_Geno${chr}.vcf.gz



  awk '{print \$1}' ${sample_GIGA} > Blood_GIGA_Common_GenoIds.txt

  vcftools --gzvcf ${GIGA} --keep Blood_GIGA_Common_GenoIds.txt  --recode  --out GIGA_Blood_Geno${chr}

  bcftools reheader -s ${sample_GIGA} -o GIGA_reheader_Blood_Geno${chr}.vcf  \
  GIGA_Blood_Geno${chr}.recode.vcf

  bgzip -f GIGA_reheader_Blood_Geno${chr}.vcf

  tabix -p vcf  GIGA_reheader_Blood_Geno${chr}.vcf.gz

 
  awk '{print 0"_"\$1}' ${sample_UAL} > Blood_UAL_Graham_Common_GenoIds.txt


  awk '{print 0"_"\$1, \$2}' ${sample_UAL} > samplelist_UAL.txt

  vcftools --gzvcf ${UAL} --keep Blood_UAL_Graham_Common_GenoIds.txt   --recode  --out UAL_Graham_Blood_Geno${chr}

  bcftools reheader -s samplelist_UAL.txt -o UAL_Graham_reheader_Blood_Geno${chr}.vcf  UAL_Graham_Blood_Geno${chr}.recode.vcf

  bgzip -f UAL_Graham_reheader_Blood_Geno${chr}.vcf

  tabix -p vcf  UAL_Graham_reheader_Blood_Geno${chr}.vcf.gz 



  awk '{print \$1}' ${sample_agvic} > Blood_ag_vic_Common_GenoIds.txt

  vcftools --gzvcf ${ag_vic} --keep Blood_ag_vic_Common_GenoIds.txt  --recode  --out ag_vic_Blood_Geno${chr}

  bcftools reheader -s ${sample_agvic} -o ag_vic_reheader_Blood_Geno${chr}.vcf  \
  ag_vic_Blood_Geno${chr}.recode.vcf

  bgzip -f ag_vic_reheader_Blood_Geno${chr}.vcf

  tabix -p vcf  ag_vic_reheader_Blood_Geno${chr}.vcf.gz



  ### --- Take care of isec -n when having the intersection of different vcf files ---

  bcftools isec -n=5 -w1  FBN_reheader_Blood_Geno${chr}.vcf.gz koch_FBN_reheader_Blood_Geno${chr}.vcf.gz GIGA_reheader_Blood_Geno${chr}.vcf.gz \
  UAL_Graham_reheader_Blood_Geno${chr}.vcf.gz ag_vic_reheader_Blood_Geno${chr}.vcf.gz \
  | bgzip -f -c > Blood_Merged_Chr${chr}.vcf.gz

  bcftools query -f '%CHROM\t%POS \n' Blood_Merged_Chr${chr}.vcf.gz > Blood_CommonVariants${chr}.txt
    

    ### -R option takes into account overlapping records. If a strict subset by position is required, add (or replace with) the -T option

  bcftools filter  --targets-file Blood_CommonVariants${chr}.txt FBN_reheader_Blood_Geno${chr}.vcf.gz \
  > FBN_Blood_CommVar_Geno${chr}.vcf

  bgzip -f FBN_Blood_CommVar_Geno${chr}.vcf

  tabix -p vcf  FBN_Blood_CommVar_Geno${chr}.vcf.gz

   
  
  bcftools filter  --targets-file Blood_CommonVariants${chr}.txt koch_FBN_reheader_Blood_Geno${chr}.vcf.gz \
  > koch_FBN_Blood_CommVar_Geno${chr}.vcf

  bgzip -f koch_FBN_Blood_CommVar_Geno${chr}.vcf

  tabix -p vcf  koch_FBN_Blood_CommVar_Geno${chr}.vcf.gz


  
  bcftools filter  --targets-file Blood_CommonVariants${chr}.txt GIGA_reheader_Blood_Geno${chr}.vcf.gz \
  > GIGA_Blood_CommVar_Geno${chr}.vcf

  bgzip -f GIGA_Blood_CommVar_Geno${chr}.vcf

  tabix -p vcf GIGA_Blood_CommVar_Geno${chr}.vcf.gz




  bcftools filter  --targets-file Blood_CommonVariants${chr}.txt UAL_Graham_reheader_Blood_Geno${chr}.vcf.gz \
  > UAL_Graham_Blood_CommVar_Geno${chr}.vcf

  bgzip -f UAL_Graham_Blood_CommVar_Geno${chr}.vcf

  tabix -p vcf UAL_Graham_Blood_CommVar_Geno${chr}.vcf.gz


   
  bcftools filter  --targets-file Blood_CommonVariants${chr}.txt ag_vic_reheader_Blood_Geno${chr}.vcf.gz\
  > ag_vic_Blood_CommVar_Geno${chr}.vcf

  bgzip -f ag_vic_Blood_CommVar_Geno${chr}.vcf

  tabix -p vcf ag_vic_Blood_CommVar_Geno${chr}.vcf.gz




  ## bcftools merge -R Blood_CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > Blood_genotypesChr${chr}.vcf

  bcftools merge *CommVar_Geno${chr}.vcf.gz > Blood_genotypesChr${chr}.vcf


  bgzip -f Blood_genotypesChr${chr}.vcf

  bcftools query -l Blood_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}



/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID
2. Extract SNPs common across all partners
3. Filter common varaints from each partner genotype data
4. # bcftools isec -n=2 -w1 extract and write records from first vcf file shared by all other vcfs using exact allele match #
*/


process tissuewise_extractGenotype_MammaryGland {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools'

 input:
      tuple val(chr), file(FBN),   file(LUKE)
      
      file(sample_FBN_HL)
      
      file(sample_FBN_VL)

      file(sample_LUKE)

   

 output:

      tuple val(chr), file ("MammaryGland_genotypesChr${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """
  ### LUKE data set has duplicate RNA seq sample for a genotype sample which has to be collected ##


  cp $sample_LUKE temp

for i in {1..17}
do
  awk '!seen[\$1]++' temp  > MammaryGland_RNA_WGS_LUKE_CorresID_set\${i}.txt
 

  awk '{print \$1}' MammaryGland_RNA_WGS_LUKE_CorresID_set\${i}.txt > MammaryGland_LUKE_Common_GenoIds.txt


  vcftools --gzvcf ${LUKE} --keep MammaryGland_LUKE_Common_GenoIds.txt  --recode  --out MG_LUKE_Common_GenoIds_set\${i}_Geno${chr}


   bcftools reheader -s MammaryGland_RNA_WGS_LUKE_CorresID_set\${i}.txt  -o LUKE_reheader_Mammary_Geno${chr}_set\${i}.vcf  \
   MG_LUKE_Common_GenoIds_set\${i}_Geno${chr}.recode.vcf

   bgzip -f LUKE_reheader_Mammary_Geno${chr}_set\${i}.vcf 

   tabix -p vcf LUKE_reheader_Mammary_Geno${chr}_set\${i}.vcf.gz


  awk 'NR==FNR {a[\$2]=\$2;next} !(\$2 in a) {print }' MammaryGland_RNA_WGS_LUKE_CorresID_set\${i}.txt temp > temp2

  mv temp2 temp

 done



  bcftools merge LUKE_reheader_Mammary_Geno${chr}_set*.vcf.gz -o LUKE_reheader_Mammary_Geno${chr}.vcf


  bgzip -f LUKE_reheader_Mammary_Geno${chr}.vcf

  tabix -p vcf  LUKE_reheader_Mammary_Geno${chr}.vcf.gz



  ### FBN data set has with the same genotype HL and VL duplicate RNA seq sample  which has to be collected ##

  
  awk '{print \$1}' ${sample_FBN_HL} > Mammary_FBN_HL_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Mammary_FBN_HL_Common_GenoIds.txt  --recode  --out FBN_HL_Mammary_Geno${chr}


  bcftools reheader -s ${sample_FBN_HL} -o FBN_HL_reheader_Mammary_Geno${chr}.vcf  \
  FBN_HL_Mammary_Geno${chr}.recode.vcf


  bgzip -f  FBN_HL_reheader_Mammary_Geno${chr}.vcf 

  tabix -p vcf FBN_HL_reheader_Mammary_Geno${chr}.vcf.gz

  awk '{print \$1}' ${sample_FBN_VL} > Mammary_FBN_VL_Common_GenoIds.txt

  vcftools --gzvcf ${FBN} --keep Mammary_FBN_VL_Common_GenoIds.txt  --recode  --out FBN_VL_Mammary_Geno${chr}


  bcftools reheader -s ${sample_FBN_VL} -o FBN_VL_reheader_Mammary_Geno${chr}.vcf  \
  FBN_VL_Mammary_Geno${chr}.recode.vcf


  bgzip -f  FBN_VL_reheader_Mammary_Geno${chr}.vcf 

  tabix -p vcf FBN_VL_reheader_Mammary_Geno${chr}.vcf.gz

  
  bcftools merge FBN_HL_reheader_Mammary_Geno${chr}.vcf.gz FBN_VL_reheader_Mammary_Geno${chr}.vcf.gz \
  -o FBN_reheader_Mammary_Geno${chr}.vcf

  bgzip -f FBN_reheader_Mammary_Geno${chr}.vcf

  tabix -p vcf  FBN_reheader_Mammary_Geno${chr}.vcf.gz


  ### --- Take care of isec -n when having the intersection of different vcf files ---

  bcftools isec -n=2 -w1  LUKE_reheader_Mammary_Geno${chr}.vcf.gz  FBN_reheader_Mammary_Geno${chr}.vcf.gz  \
   | bgzip -f -c > MammaryGland_Merged_Chr${chr}.vcf.gz

  bcftools query -f '%CHROM\t%POS \n' MammaryGland_Merged_Chr${chr}.vcf.gz > MammaryGland_CommonVariants${chr}.txt
    

    ### -R option takes into account overlapping records. If a strict subset by position is required, add (or replace with) the -T option

  bcftools filter  --targets-file MammaryGland_CommonVariants${chr}.txt FBN_reheader_Mammary_Geno${chr}.vcf.gz \
  > FBN_Mammary_CommVar_Geno${chr}.vcf

  bgzip -f FBN_Mammary_CommVar_Geno${chr}.vcf

  tabix -p vcf  FBN_Mammary_CommVar_Geno${chr}.vcf.gz

   
  bcftools filter  --targets-file MammaryGland_CommonVariants${chr}.txt LUKE_reheader_Mammary_Geno${chr}.vcf.gz \
  > LUKE_Mammary_CommVar_Geno${chr}.vcf

  bgzip -f LUKE_Mammary_CommVar_Geno${chr}.vcf

  tabix -p vcf  LUKE_Mammary_CommVar_Geno${chr}.vcf.gz


  ## bcftools merge -R MammaryGland_CommonVariants${chr}.txt  *CommVar_Geno${chr}.vcf.gz > MammaryGland_genotypesChr${chr}.vcf

  bcftools merge *CommVar_Geno${chr}.vcf.gz > MammaryGland_genotypesChr${chr}.vcf


  bgzip -f MammaryGland_genotypesChr${chr}.vcf

  bcftools query -l MammaryGland_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}


/*

1. Extract genotype data from each individual partner and rename the samples based on RNAseq samples ID


*/


process tissuewise_extractGenotype_Milk {
 tag " on chromosome $chr"
 publishDir "${params.outputGeno}/eQTLGenoData", mode:'copy'
 container 'praveen/qtltools'

 input:
      tuple val(chr),  file(ag_vic)

      file(sample_agvic)

 output:

      tuple val(chr), file ("Milk_genotypesChr${chr}.vcf.gz")

      file("GenoSamplesChr${chr}.txt")
      

 script:

 """

  awk '{print \$1}' ${sample_agvic} > Milk_ag_vic_Common_GenoIds.txt

  vcftools --gzvcf ${ag_vic} --keep Milk_ag_vic_Common_GenoIds.txt  --recode  --out ag_vic_Milk_Geno${chr}

  bcftools reheader -s ${sample_agvic} -o Milk_genotypesChr${chr}.vcf  \
  ag_vic_Milk_Geno${chr}.recode.vcf


  bgzip -f Milk_genotypesChr${chr}.vcf

  bcftools query -l Milk_genotypesChr${chr}.vcf.gz | awk '{print \$1".junc"}' > GenoSamplesChr${chr}.txt

 """

}