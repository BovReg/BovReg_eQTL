
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

process transcriptCountsPCA_trans {
 tag "on chromosome ${chr}"
 publishDir "${params.outdir}", mode:'copy'
 container 'praveenchitneedi/qtltoolkit:v1.0.1'

 input:
  file (phenotype) 

  tuple val(chr),file (genotype) ,file (genoCovfile)

  val(phenotype_PCs_trans)


 output:

 tuple val(chr),file ("Geno10_Chr${chr}_PhenoAllChr_transcript.pca*") , emit:  cis_covariates_transcript_ch

 tuple val(chr), file ("QTLphenoChr${chr}_transcript_Allchr.bed.*") , emit: phenotype_transcript_Bed_ch

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


 ## Extract bed file with all chromosomes to perform trans-eQTL across all the chromosomal phenotypes ##

 sed  's/Chr/#CHROM/g' transcript_count_Matrices_CommonSamples.tsv > QTLphenoChr${chr}_transcript_Allchr.bed

 ## awk: line 1: runaway string constant " ...  should use \\n in your original string.


awk -F'\t' '{for(i=1;i<=6;i++){printf "%s ", \$i}; printf "\\n"}' QTLphenoChr${chr}_transcript_Allchr.bed > bed_pos

awk -F'\t' '{for(i=7;i<=NF;i++){printf "%s ", \$i}; printf "\\n"}' QTLphenoChr${chr}_transcript_Allchr.bed > bed_data
 
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


 paste -d' ' bed_pos temp | awk -v OFS="\t" '\$1=\$1'  > QTLphenoChr${chr}_transcript_Allchr.bed

 awk 'FNR==1{print}' QTLphenoChr${chr}_transcript_Allchr.bed > temp1

 awk 'FNR>1{print}' QTLphenoChr${chr}_transcript_Allchr.bed | awk -F "\t" '{ if(\$1 >= 1 && \$1 <= 29){ print } }' | sort -n -k1,1 -k2,2 > temp2

 cat temp1 temp2 > QTLphenoChr${chr}_transcript_Allchr.bed
 
 rm bed_pos temp temp2 bed_data

 bgzip QTLphenoChr${chr}_transcript_Allchr.bed && tabix -p bed QTLphenoChr${chr}_transcript_Allchr.bed.gz

 QTLtools pca --bed QTLphenoChr${chr}_transcript_Allchr.bed.gz --scale --center --out pheno_transcript_QTL.Allchr

 ## NOTE: The default PCs is 10, but user can modify it by changing 'NR<11{ print }' in json file 
 awk 'NR<=${phenotype_PCs_trans}+1{ print }' pheno_transcript_QTL.Allchr.pca > pheno_transcript_QTL.Allchr_top.pca

  ### Extract common samples and merge common columns from phenotype and genotype covariates

 awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} } ' pheno_transcript_QTL.Allchr_top.pca ${genoCovfile} \
    | awk -v OFS="\t" '\$1=\$1' > Geno10_Chr${chr}_PhenoAllChr_transcript.pca

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
  }' Geno10_Chr${chr}_PhenoAllChr_transcript.pca > xx

  ## Re-order columns according to vcf samples

   gawk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx sample_listpca.txt  > xy 

 ## Re-transpose file and remove genotype covariates 
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
  }' xy | sed '7,11d;17,21d'  > Geno10_Chr${chr}_PhenoAllChr_transcript.pca

 bgzip Geno10_Chr${chr}_PhenoAllChr_transcript.pca 
  
 gawk  'NR==1{for (i=2;i<=NF;i++) print \$i}'  pheno_transcript_QTL.Allchr.pca > samplefromRNAseq.txt

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

process geneCountsPCA_trans {
 tag "on chromosome ${chr}"
 publishDir "${params.outdir}", mode:'copy'
 container 'praveenchitneedi/qtltoolkit:v1.0.1'

 input:
 file (phenotype) 

 tuple val(chr),  file (genotype), file (genoCovfile)

 val(phenotype_PCs_trans)

 output:

 
 tuple val(chr), path("phenoChr${chr}_gene_QTL_Allchr.bed.*") , emit: phenotype_gene_Bed_ch

 tuple val(chr), file ("Geno10_Chr${chr}_PhenoAllChr_Gene.pca*") , emit:  cis_covariates_gene_ch

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

## Extract bed file with all chromosomes to perform trans-eQTL across all the chromosomal phenotypes ##

 sed  's/Chr/#CHROM/g' Gene_count_Matrices_CommonSamples.tsv >  phenoChr${chr}_gene_QTL_Allchr.bed

awk -F'\t' '{for(i=1;i<=6;i++){printf "%s ", \$i}; printf "\\n"}' phenoChr${chr}_gene_QTL_Allchr.bed > bed_pos

 awk -F'\t' '{for(i=7;i<=NF;i++){printf "%s ", \$i}; printf "\\n"}' phenoChr${chr}_gene_QTL_Allchr.bed > bed_data 
 

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

 paste -d' ' bed_pos temp | awk -v OFS="\t" '\$1=\$1'  > phenoChr${chr}_gene_QTL_Allchr.bed
 
 awk 'FNR==1{print}' phenoChr${chr}_gene_QTL_Allchr.bed > temp1

 awk 'FNR>1{print}' phenoChr${chr}_gene_QTL_Allchr.bed | awk -F "\t" '{ if(\$1 >= 1 && \$1 <= 29){ print } }' | sort -n -k1,1 -k2,2 > temp2

 cat temp1 temp2 > phenoChr${chr}_gene_QTL_Allchr.bed

 rm bed_pos temp bed_data temp2 

 bgzip phenoChr${chr}_gene_QTL_Allchr.bed && tabix -p bed phenoChr${chr}_gene_QTL_Allchr.bed.gz

 QTLtools pca --bed phenoChr${chr}_gene_QTL_Allchr.bed.gz --scale --center --out pheno_gene_QTL_chr

 ## NOTE: The default PCs is 10, but user can modify it by changing 'NR<11{ print }' in json file 
 awk 'NR<=${phenotype_PCs_trans}+1{ print }' pheno_gene_QTL_chr.pca > pheno_gene_QTL_Allchr_top.pca

 ### Extract common samples and merge common columns from phenotype and genotype covariates

awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} }' pheno_gene_QTL_Allchr_top.pca ${genoCovfile} \
   | awk -v OFS="\t" '\$1=\$1' > Geno10_Chr${chr}_PhenoAllChr_Gene.pca

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
  }' Geno10_Chr${chr}_PhenoAllChr_Gene.pca > xx2

  ## Re-order columns according to vcf samples

   gawk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx2 sample_listpca.txt  > xy2 

 ## Re-transpose file and remove genotype covariates 
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
  }' xy2 | sed '7,11d;17,21d'  > Geno10_Chr${chr}_PhenoAllChr_Gene.pca


 bgzip Geno10_Chr${chr}_PhenoAllChr_Gene.pca 

 awk  'NR==1{for (i=2;i<=NF;i++) print \$i}' pheno_gene_QTL_chr.pca > samplefromRNAseq.txt


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

process rnaSplicePCS_trans {
  tag "on chromosome ${chr}"
  publishDir "${params.outdir}", mode:'copy'
  container 'praveenchitneedi/qtltoolkit:v1.0.1'
   input:

   file (covariate_s) 

   file (phenotype_s)

   tuple val(chr),  file (genotype),  file (genoCovfile)

   output:

   tuple val(chr), path ("Geno10_Pheno10_Splice_Chr${chr}.pca.*") , emit: cis_covariates_splice_ch

   tuple val(chr), path ("phenochr${chr}_splice_QTL_Allchr.bed*") , emit: phenotype_splice_Bed_ch

   tuple val(chr), path  ("Geno_Chr${chr}_splice.vcf.*") , emit: genotypeQTL_splice_ch

  script:

//## Merge the data of all the phenotype to perform trans-eQTL with all the phenotypes acoss genome for each variant ##
// for phenosplice in ${phenotype_s}
//  do

 //  cat \$phenosplice | awk 'FNR>1{print}' >> pheno_splice_QTL_chr.bed 

  // cat \$phenosplice | awk 'FNR==1{print}' > pheno_splice_QTL_chr_header.bed
 //done

//sort -n -k1,1 -k2,2 pheno_splice_QTL_chr.bed > pheno_splice_QTL_chr_sorted.bed

//cat pheno_splice_QTL_chr_header.bed pheno_splice_QTL_chr_sorted.bed > phenochr${chr}_splice_QTL_Allchr.bed



 """

 bcftools query -l ${genotype} > samplelist

 sed '1iid' samplelist > sample_listpca.txt

 awk ' NR == FNR { header[\$0]=1; next } 
    FNR == 1 { 
    for (i=1; i<=NF; i++) if (\$i in header) wanted[i]=1 }
    { for (i=1; i<=NF; i++) if (i in wanted) printf "%s ", \$i;  print "" }' samplelist ${phenotype_s} > xx

 awk '{for(i=1;i<=6;i++) printf \$i" "; print ""}' ${phenotype_s} > xy

 paste -d ' ' xy xx | awk -v OFS="\t" '\$1=\$1' | sed 's/Chr/#CHROM/g'  >  phenochr${chr}_splice_QTL_Allchr.bed

 awk -F'\t' '{for(i=1;i<=6;i++){printf "%s ", \$i}; printf "\\n"}' phenochr${chr}_splice_QTL_Allchr.bed > bed_pos

 awk -F'\t' '{for(i=7;i<=NF;i++){printf "%s ", \$i}; printf "\\n"}' phenochr${chr}_splice_QTL_Allchr.bed > bed_data 
 
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

 paste -d' ' bed_pos temp | awk -v OFS="\t" '\$1=\$1'  > phenochr${chr}_splice_QTL_Allchr.bed
 
 rm bed_pos temp bed_data


 bgzip  phenochr${chr}_splice_QTL_Allchr.bed && tabix -p bed phenochr${chr}_splice_QTL_Allchr.bed.gz

 sed 's/comm_geno_//g'  $covariate_s | awk 'NR==1;NR>1{\$1="splicePC_"\$1 ;print}' | awk -v OFS=' ' '{\$1=\$1}1' >  ${covariate_s}_MOD.pca

 sed  's/SampleID/id/g' ${genoCovfile} > genCovfile_Mod.txt

  ### Extract common samples and merge common columns from phenotype and genotype covariates

 awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} } '  ${covariate_s}_MOD.pca genCovfile_Mod.txt \
    > Geno10_Pheno10_Splice_Chr${chr}.pca

  awk 'NR==1{for (i=2; i<=NF; i++) print \$i}' ${covariate_s}_MOD.pca > Ind_leafcutter.txt

 ## Transpose file and remove genotype covariates for testing
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
  }' Geno10_Pheno10_Splice_Chr${chr}.pca > xx1

  ## Re-order columns according to vcf samples

   gawk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx1 sample_listpca.txt  > xy2 

 ## Re-transpose file and remove genotype covariates
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
  }' xy2  | sed '7,11d;17,21d' > Geno10_Pheno10_Splice_Chr${chr}.pca



  bgzip Geno10_Pheno10_Splice_Chr${chr}.pca 

  ## zcat $genotype  | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > temp

  vcftools --gzvcf $genotype --keep Ind_leafcutter.txt --maf 0.05  --recode --out Geno_Chr${chr}_splice  

  mv Geno_Chr${chr}_splice.recode.vcf Geno_Chr${chr}_splice.vcf
   
  bgzip Geno_Chr${chr}_splice.vcf && tabix -p vcf Geno_Chr${chr}_splice.vcf.gz

   """
}