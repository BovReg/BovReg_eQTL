
/*
process TranscriptCountsPCA: The continer used for this process is preperared with the dockerfile
/home/chitneedi/Disk2chitneedi/NextFlow/QTLtools_Dockerfile/Dockerfile
 This dockerfile has the installations for software QTLtools, fatqtl, tabix, vcftools, samtools, bcftools and htslib
 from biocontainers "https://github.com/BioContainers/containers"

 - extract header and substring the sample id
 - extract the non-header data and split the phenotype id , emit: Chr, start, end, pid,gid, strand and samples
 - merge the header and data , emit: .bed file
 - compress the .bed file , emit: bgzip and create tabix for the compressed file
 - perform pca analysis with the compressed .bed file uing QTLtools
 - Extract the top 10 phentoype PCs.
 - Concatenate the top 10 phenotype and top 10 genotype PCs
 - Extract samples common to phenotype data from genotype vcf file using vcftools.
 - create a new variable for each chromosome excluding the .vcf extension.


*/

process transcriptCountsPCA {
 publishDir "${params.outputQTL}", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
  file (phenotype) 

  tuple val(chr),file (genotype) ,file (genoCovfile)

  val(phenotype_PCs_cis)


 output:

 tuple val(chr),file ("Geno10_Pheno10_transcript_Chr${chr}.pca*") , emit:  cis_covariates_transcript_ch

 tuple val(chr),file ("QTLpheno_transcript_chr${chr}.bed.*") , emit: phenotype_transcript_Bed_ch

 tuple val(chr), file ("Geno_Chr${chr}_transcript.vcf.*"), emit: genotypeQTL_transcript_ch

 script:

 """

  bcftools query -l ${genotype} > samplelist

  sed '1iSampleID' samplelist > sample_listpca.txt


  awk ' NR == FNR { header[\$0]=1; next } 
  FNR == 1 { 
    for (i=1; i<=NF; i++) if (\$i in header) wanted[i]=1 }
    { for (i=1; i<=NF; i++) if (i in wanted) printf "%s ", \$i;  print "" }' samplelist ${phenotype} > xx

  awk '{for(i=1;i<=6;i++) printf \$i" "; print ""}' ${phenotype} > xy

  paste -d ' ' xy xx | awk -v OFS="\t" '\$1=\$1' >  transcript_count_Matrices_CommonSamples.tsv

 ## Extract bed file for each chromosome to perform chromsomal wise cis-eQTL analysis ##
 
 awk -v OFS='\t' 'NR==1 {\$1="#Chr";print }' transcript_count_Matrices_CommonSamples.tsv  > xy

 awk -v OFS='\t' 'NR>1 && \$1==${chr} {print}' transcript_count_Matrices_CommonSamples.tsv | sort -k1,1d -k2,2n -k3,3n > xx

 cat xy xx > QTLpheno_transcript_chr${chr}.bed

 sed -i 's/#Chr/#CHROM/g' QTLpheno_transcript_chr${chr}.bed

 ## awk: line 1: runaway string constant " ...  should use \\n in your original string.

 awk -F'\t' '{for(i=1;i<=6;i++){printf "%s ", \$i}; printf "\\n"}' QTLpheno_transcript_chr${chr}.bed > bed_pos

 awk -F'\t' '{for(i=7;i<=NF;i++){printf "%s ", \$i}; printf "\\n"}' QTLpheno_transcript_chr${chr}.bed > bed_data 
 
 ## Transpose file 
  awk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' bed_data > xx


 awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx samplelist  > xy 

  ## Re-transpose file
 awk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' xy > temp

 paste -d' ' bed_pos temp | awk -v OFS="\t" '\$1=\$1'  > QTLpheno_transcript_chr${chr}.bed
 
 rm bed_pos temp bed_data

 bgzip QTLpheno_transcript_chr${chr}.bed && tabix -p bed QTLpheno_transcript_chr${chr}.bed.gz

 QTLtools pca --bed QTLpheno_transcript_chr${chr}.bed.gz --scale --center --out pheno_transcript_QTL.chr${chr}

 ## NOTE: The default PCs is 10, but user can modify it by changing 'NR<11{ print }' in json file 
 awk 'NR<=${phenotype_PCs_cis}+1{ print }' pheno_transcript_QTL.chr${chr}.pca > pheno_transcript_QTL.chr${chr}_top.pca

  ### Extract common samples and merge common columns from phenotype and genotype covariates

 awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} } ' pheno_transcript_QTL.chr${chr}_top.pca ${genoCovfile} \
    | awk -v OFS="\t" '\$1=\$1' > Geno10_Pheno10_transcript_Chr${chr}.pca

  ## Transpose file 
  gawk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' Geno10_Pheno10_transcript_Chr${chr}.pca > xx

  ## Re-order columns according to vcf samples

   gawk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx sample_listpca.txt  > xy 


 ## Re-transpose file
 gawk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' xy > Geno10_Pheno10_transcript_Chr${chr}.pca

 bgzip Geno10_Pheno10_transcript_Chr${chr}.pca 

 gawk  'NR==1{for (i=2;i<=NF;i++) print \$i}' pheno_transcript_QTL.chr${chr}.pca > samplefromRNAseq.txt

 ##zcat $genotype | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.1/'  > temp

 vcftools --gzvcf $genotype --keep samplefromRNAseq.txt --maf 0.05 --recode  --out Geno_Chr${chr}_transcript
 
 mv Geno_Chr${chr}_transcript.recode.vcf Geno_Chr${chr}_transcript.vcf

 bgzip Geno_Chr${chr}_transcript.vcf && tabix -p vcf Geno_Chr${chr}_transcript.vcf.gz

 rm  xx xy

"""
}




/*
process GeneCountsPCA: The continer used for this process is preperared with the dockerfile
/home/chitneedi/Disk2chitneedi/NextFlow/QTLtools_Dockerfile/Dockerfile
 This dockerfile has the installations for software QTLtools, fatqtl, tabix, samtools, vcftools, bcftools and htslib
 from biocontainers "https://github.com/BioContainers/containers"

 - extract header and substring the sample id
 - extract the non-header data and split the phenotype id , emit: Chr, start, end, pid,gid, strand and samples
 - merge the header and data , emit: .bed file
 - compress the .bed file , emit: bgzip and create tabix for the compressed file
 - perform pca analysis with the compressed .bed file uing QTLtools
 - Extract the top 10 phentoype PCs.
 - Concatenate the top 10 phenotype and top 10 genotype PCs
 - Extract samples common to phenotype data from genotype vcf file using vcftools.
 - create a new variable for each chromosome excluding the .vcf extension.

*/

process geneCountsPCA {
 publishDir "${params.outputQTL}", mode:'copy'
 container 'praveen/qtltools1.3'

 input:
 file (phenotype) 

 tuple val(chr),  file (genotype), file (genoCovfile)

 val(phenotype_PCs_cis)

 output:

 
 tuple val(chr), path("pheno_gene_QTL_chr${chr}.bed.*") , emit: phenotype_gene_Bed_ch

 tuple val(chr), file ("Geno10_Pheno10_Geno_Chr${chr}.pca*") , emit:  cis_covariates_gene_ch

 tuple val(chr), path ("Geno_Chr${chr}_gene.vcf.*"),  emit: genotypeQTL_gene_ch 

 
 script:

  """

  bcftools query -l ${genotype} > samplelist

  sed '1iSampleID' samplelist > sample_listpca.txt

  awk ' NR == FNR { header[\$0]=1; next } 
  FNR == 1 { 
    for (i=1; i<=NF; i++) if (\$i in header) wanted[i]=1 }
    { for (i=1; i<=NF; i++) if (i in wanted) printf "%s ", \$i;  print "" }' samplelist ${phenotype}  > xx

  awk '{for(i=1;i<=6;i++) printf \$i" "; print ""}' ${phenotype} > xy

  paste -d ' ' xy xx | awk -v OFS="\t" '\$1=\$1' >  Gene_count_Matrices_CommonSamples.tsv


 awk -v OFS='\t' 'NR==1 {\$1="#Chr";print }' Gene_count_Matrices_CommonSamples.tsv  > xy

 ## Extract bed file for each chromosome to perform chromsomal wise cis-eQTL analysis ##

 awk -v OFS='\t' 'NR>1 && \$1==${chr} {print}'  Gene_count_Matrices_CommonSamples.tsv | sort -k1,1d -k2,2n -k3,3n > xx

 cat xy xx > pheno_gene_QTL_chr${chr}.bed

 sed -i 's/#Chr/#CHROM/g' pheno_gene_QTL_chr${chr}.bed

awk -F'\t' '{for(i=1;i<=6;i++){printf "%s ", \$i}; printf "\\n"}' pheno_gene_QTL_chr${chr}.bed > bed_pos

 awk -F'\t' '{for(i=7;i<=NF;i++){printf "%s ", \$i}; printf "\\n"}' pheno_gene_QTL_chr${chr}.bed > bed_data 
 
 ## Transpose file 
  awk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' bed_data > xx


 awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx samplelist  > xy 

  ## Re-transpose file
 awk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' xy > temp

 paste -d' ' bed_pos temp | awk -v OFS="\t" '\$1=\$1'  > pheno_gene_QTL_chr${chr}.bed
 
 rm bed_pos temp bed_data


 bgzip pheno_gene_QTL_chr${chr}.bed && tabix -p bed pheno_gene_QTL_chr${chr}.bed.gz

 QTLtools pca --bed pheno_gene_QTL_chr${chr}.bed.gz --scale --center --out pheno_gene_QTL_chr${chr}

 ## NOTE: The default PCs is 10, but user can modify it by changing 'NR<11{ print }' in json file
 awk 'NR<${phenotype_PCs_cis}+1{ print }' pheno_gene_QTL_chr${chr}.pca > pheno_gene_QTL_chr${chr}_top.pca

 ### Extract common samples and merge common columns from phenotype and genotype covariates

awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} }' pheno_gene_QTL_chr${chr}_top.pca ${genoCovfile} \
   | awk -v OFS="\t" '\$1=\$1' > Geno10_Pheno10_Geno_Chr${chr}.pca

## Transpose file 
  gawk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' Geno10_Pheno10_Geno_Chr${chr}.pca > xx2

  ## Re-order columns according to vcf samples

   gawk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx2 sample_listpca.txt  > xy2 

 ## Re-transpose file
 gawk '
  { 
      for (i=1; i<=NF; i++)  {
          a[NR,i] = \$i
      }
  }
  NF>p { p = NF }
  END {    
      for(j=1; j<=p; j++) {
          str=a[1,j]
          for(i=2; i<=NR; i++){
              str=str" "a[i,j];
          }
          print str
      }
  }' xy2 > Geno10_Pheno10_Geno_Chr${chr}.pca


 bgzip Geno10_Pheno10_Geno_Chr${chr}.pca 

 awk  'NR==1{for (i=2;i<=NF;i++) print \$i}' pheno_gene_QTL_chr${chr}.pca > samplefromRNAseq.txt


 ## zcat $genotype | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.1/'  > temp
 
 vcftools --gzvcf $genotype  --keep samplefromRNAseq.txt  --maf 0.05  --recode --out Geno_Chr${chr}_gene

 mv Geno_Chr${chr}_gene.recode.vcf Geno_Chr${chr}_gene.vcf

 bgzip Geno_Chr${chr}_gene.vcf && tabix -p vcf Geno_Chr${chr}_gene.vcf.gz

 rm  xx xy

"""
}



/*

 process RNAsplicePCS: The continer used for this process is preperared with the dockerfile
 /home/chitneedi/Disk2chitneedi/NextFlow/QTLtools_Dockerfile/Dockerfile
 This dockerfile has the installations for software QTLtools, fatqtl, tabix, vcftools, samtools, bcftools and htslib
 from biocontainers "https://github.com/BioContainers/containers"
 
 - The splice junctions phenotype data was used a input which was created using "process leafcutter_cluster_junctions in eQTL_RNAseq_MappingQunatif_01.nf"
 - extract header and substring the sample id
 - extract the non-header data and split the phenotype id , emit: Chr, start, end, pid,gid, strand and samples
 - merge the header and data , emit: .bed file
 - compress the .bed file , emit: bgzip and create tabix for the compressed file
 - Concatenate the phenotype and top 10 genotype PCs
 - Extract samples common to phenotype data from genotype vcf file using vcftools.
 - create a new variable for each chromosome excluding the .vcf extension.


*/

process rnaSplicePCS {

  publishDir "${params.outputQTL}", mode:'copy'
  container 'praveen/qtltools1.3'
   input:

   file (covariate_s) 

   tuple val(chr), file (phenotype_s),  file (genotype),  file (genoCovfile)

   output:

   tuple val(chr), path ("Geno10_Pheno10_Splice_Chr${chr}.pca.*") , emit: cis_covariates_splice_ch

   tuple val(chr), path ("${phenotype_s}.bed.*") , emit: phenotype_splice_Bed_ch

   tuple val(chr), path  ("Geno_Chr${chr}_splice.vcf.*") , emit: genotypeQTL_splice_ch

  script:

 """

    bcftools query -l ${genotype} > samplelist

    sed '1iid' samplelist > sample_listpca.txt


     sed 's/comm_//g' $phenotype_s | awk -F'\t' 'NR>1{print \$1,\$2,\$3,\$4}' | awk -F'_' '{print \$1"_"\$2"\t"\$3}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$1":"\$2"-"\$3,\$4,\$5}' | sed '1i#Chr\tstart\tend\tpid\tgid\tstrand' > ${phenotype_s}_MOD.bed.header

     sed 's/comm_//g' $phenotype_s  | awk  '{for(i=5;i<=NF;i++) printf \$i"\t"; print ""}'  > ${phenotype_s}_MOD.bed.counts

     paste -d '\t' ${phenotype_s}_MOD.bed.header ${phenotype_s}_MOD.bed.counts | awk -v OFS='\t' 'NR==1 {print}' > xx

     paste -d'\t' ${phenotype_s}_MOD.bed.header ${phenotype_s}_MOD.bed.counts |  awk -v OFS='\t' 'NR>1 {print}' | sort -k1,1d -k2,2n -k3,3n > xy

     cat xx xy | sed 's/geno_//g' > ${phenotype_s}.bed

     rm ${phenotype_s}_MOD.bed.header ${phenotype_s}_MOD.bed.counts

     awk -F'\t' '{for(i=1;i<=6;i++){printf "%s ", \$i}; printf "\\n"}' ${phenotype_s}.bed > bed_pos

     awk -F'\t' '{for(i=7;i<=NF;i++){printf "%s ", \$i}; printf "\\n"}' ${phenotype_s}.bed > bed_data 
     
     ## Transpose file 
      awk '
      { 
          for (i=1; i<=NF; i++)  {
              a[NR,i] = \$i
          }
      }
      NF>p { p = NF }
      END {    
          for(j=1; j<=p; j++) {
              str=a[1,j]
              for(i=2; i<=NR; i++){
                  str=str" "a[i,j];
              }
              print str
          }
      }' bed_data > xx


     awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx samplelist  > xy 

      ## Re-transpose file
     awk '
      { 
          for (i=1; i<=NF; i++)  {
              a[NR,i] = \$i
          }
      }
      NF>p { p = NF }
      END {    
          for(j=1; j<=p; j++) {
              str=a[1,j]
              for(i=2; i<=NR; i++){
                  str=str" "a[i,j];
              }
              print str
          }
      }' xy > temp

     paste -d' ' bed_pos temp | awk -v OFS="\t" '\$1=\$1'  > ${phenotype_s}.bed
     
     rm bed_pos temp bed_data


     bgzip  ${phenotype_s}.bed && tabix -p bed ${phenotype_s}.bed.gz

     sed 's/comm_geno_//g'  $covariate_s | awk 'NR==1;NR>1{\$1="splicePC_"\$1 ;print}' | awk -v OFS=' ' '{\$1=\$1}1' >  ${covariate_s}_MOD.pca

     sed  's/SampleID/id/g' ${genoCovfile} > genCovfile_Mod.txt

      ### Extract common samples and merge common columns from phenotype and genotype covariates

     awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
     FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
      FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} } '  ${covariate_s}_MOD.pca genCovfile_Mod.txt \
        > Geno10_Pheno10_Splice_Chr${chr}.pca

      awk 'NR==1{for (i=2; i<=NF; i++) print \$i}' ${covariate_s}_MOD.pca > Ind_leafcutter.txt

     ## Transpose file 
      gawk '
      { 
          for (i=1; i<=NF; i++)  {
              a[NR,i] = \$i
          }
      }
      NF>p { p = NF }
      END {    
          for(j=1; j<=p; j++) {
              str=a[1,j]
              for(i=2; i<=NR; i++){
                  str=str" "a[i,j];
              }
              print str
          }
      }' Geno10_Pheno10_Splice_Chr${chr}.pca > xx

      ## Re-order columns according to vcf samples

       gawk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx sample_listpca.txt  > xy 

     ## Re-transpose file
     gawk '
      { 
          for (i=1; i<=NF; i++)  {
              a[NR,i] = \$i
          }
      }
      NF>p { p = NF }
      END {    
          for(j=1; j<=p; j++) {
              str=a[1,j]
              for(i=2; i<=NR; i++){
                  str=str" "a[i,j];
              }
              print str
          }
      }' xy > Geno10_Pheno10_Splice_Chr${chr}.pca


      bgzip Geno10_Pheno10_Splice_Chr${chr}.pca 

      ## zcat $genotype  | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > temp

      vcftools --gzvcf $genotype --keep Ind_leafcutter.txt --maf 0.05  --recode --out Geno_Chr${chr}_splice  

      mv Geno_Chr${chr}_splice.recode.vcf Geno_Chr${chr}_splice.vcf
       
      bgzip Geno_Chr${chr}_splice.vcf && tabix -p vcf Geno_Chr${chr}_splice.vcf.gz

      rm  xx xy

   """
}
