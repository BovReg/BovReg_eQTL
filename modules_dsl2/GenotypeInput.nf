

/*
process extractCommmonSamples:
       -This process extract the samples from genotpye common to RNAseq
       - The ids are changed  to the RNAseq sample Ids for QTL mapping
      
*/

process extractCommmonSamples {

  container 'praveen/qtltools'

  input:

      tuple val (chr), file(genotype)

      file (sampleinfo_FBN_RNA)

      file (liversamples_FBN)

  output:

    tuple val(chr), file ("Chr${chr}_Genotype_InputQTL.vcf")

  script:

  outgeno = genotype[0].toString() - ~/(\.vcf)?(\.gz)?$/

  outgeno2 = outgeno - ~/ImpuHDktoWGS_segfam_ChromMas_NUDA_Ars1.2_NdDidierfile_1043IND_/
  """
   awk -F'\t' 'NR>1{print "0_"\$5,\$2}' $sampleinfo_FBN_RNA > WGS_FBN_sample.txt
   
   awk 'NR==1 {print}' $liversamples_FBN | tr '\t' '\n' | awk 'NR>6{print }' > Liver_samples.txt

   awk 'NR==FNR {a[\$1]=\$1;next} \$2 in a {print }' Liver_samples.txt WGS_FBN_sample.txt > Liver_WGS_weitje_sample.txt

   awk 'NR==FNR {a[\$2]=\$0;next} !(\$1 in a) {print "0_"\$1,\$1}' Liver_WGS_weitje_sample.txt Liver_samples.txt \
    | awk '{ if( length(\$1)>=6) gsub(/H0/,"H",\$1);  print  }' > Liver_WGS_ChronMas_sample.txt

   cat Liver_WGS_weitje_sample.txt Liver_WGS_ChronMas_sample.txt | awk -v OFS='\t' '{print }' > Liver_WGS_weitje_ChronMas_sample.txt 

   awk '{print \$1}' Liver_WGS_weitje_ChronMas_sample.txt > WGS_weitje_ChronMas_sampleIds.txt


   vcftools --keep WGS_weitje_ChronMas_sampleIds.txt \
   --gzvcf $genotype --recode --out Chr${chr}_Genotype

   
   bcftools query -l Chr${chr}_Genotype.recode.vcf > sampleorder.txt
 
   awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}'  Liver_WGS_weitje_ChronMas_sample.txt sampleorder.txt \
    > Liver_WGS_weitje_ChronMas_sample_reorder.txt

   bcftools reheader -s Liver_WGS_weitje_ChronMas_sample_reorder.txt \
   Chr${chr}_Genotype.recode.vcf -o Chr${chr}_Genotype_InputQTL.vcf

"""
}