
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
 publishDir "${params.outputQTL}/TestPCAs", mode:'copy'
 container 'praveen/qtltools'

 input:
  file (phenotype) 

  tuple val(chr),file (genotype) ,file (genoCovfile)


 output:

 tuple val(chr),file ("Geno10_Pheno10_transcript_Chr${chr}.pca*") , emit:  cis_covariates_transcript_ch

 tuple val(chr),file ("QTLpheno_transcript_chr${chr}.bed.*") , emit: phenotype_transcript_Bed_ch

 tuple val(chr), file ("Geno_Chr${chr}_transcript.recode.vcf.*"), emit: genotypeQTL_transcript_ch

 script:
  //outgeno = genotype[0].toString() - ~/()?(_InputQTL)?(\.vcf)?$/
  //chr = outgeno - ~/ImpuHDktoWGS_RNA_liver_Chr_/
  /*
   awk -v OFS='\t' 'NR==1 {\$1="#Chr";print }' $phenotype | awk -F'\t' -v OFS='\t' ' {for (i=1;i<=NF;i++) if (i>6) \$i=substr(\$i,1,4)  ;else \$i=\$i } END{print \$0}' > xy
   Extract the common samples between pehotype and genotype covariate file ####
   awk 'NR==1{for (i=1;i<=NF;i++) print \$i}' pheno_transcript_QTL.Allchr.pca > sampleforGenoCov.txt
   awk 'NR==FNR{arr[\$1]++;next}{for(i=1; i<=NF; i++) if (\$i in arr){a[i]++;}}  { for (i in a) printf "%s ", \$i; printf "\n"}'  sampleforGenoCov.txt  ${genoCovfile} > genoCovfile_CommonIND.txt 
   Transpose Cov file and re arrange the samples according the phenotype data and re-transpose to original form #
   awk '{ for (i = 1; i <= NF; i++) { s[i] = s[i]?s[i] FS \$i:\$i } } END { for (i in s) { print s[i] } }' genoCovfile_CommonIND.txt  > genoCovfile_CommonIND_T.txt
   awk 'NR==FNR{o[FNR]=\$1; next} {t[\$1]=\$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}'  sampleforGeno.txt genoCovfile_CommonIND_T.txt > temp
   awk '{ for (i = 1; i <= NF; i++) { s[i] = s[i]?s[i] FS \$i:\$i } } END { for (i in s) { print s[i] } }' temp | sort -nk1 > genoCovfile_CommOrdIND.txt
 */
 """

  bcftools query -l ${genotype} > samplelist

  sed '1iSampleID' samplelist > sample_listpca.txt

  awk ' NR == FNR { header[\$0]=1; next } 
  FNR == 1 { 
    for (i=1; i<=NF; i++) if (\$i in header) wanted[i]=1 }
    { for (i=1; i<=NF; i++) if (i in wanted) printf "%s ", \$i;  print "" }' samplelist ${phenotype} > xx

  awk '{for(i=1;i<=6;i++) printf \$i" "; print ""}' ${phenotype} > xy

  paste -d ' ' xy xx | awk -v OFS="\t" '\$1=\$1' >  transcript_count_Matrices_CommonSamples.tsv


 awk -v OFS='\t' 'NR==1 {\$1="#Chr";print }' transcript_count_Matrices_CommonSamples.tsv  > xy

 awk -v OFS='\t' 'NR>1 && \$1==${chr} {print}' transcript_count_Matrices_CommonSamples.tsv | sort -k1,1d -k2,2n -k3,3n > xx

 cat xy xx > QTLpheno_transcript_chr${chr}.bed

 sed -i 's/#Chr/#CHROM/g' QTLpheno_transcript_chr${chr}.bed

 bgzip QTLpheno_transcript_chr${chr}.bed && tabix -p bed QTLpheno_transcript_chr${chr}.bed.gz

 QTLtools pca --bed QTLpheno_transcript_chr${chr}.bed.gz --scale --center --out pheno_transcript_QTL.chr${chr}

 awk 'NR<12{ print }' pheno_transcript_QTL.chr${chr}.pca > pheno_transcript_QTL.chr${chr}_top.pca


  ### Extract common samples and merge common columns from phenotype and genotype covariates

 awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} } ' pheno_transcript_QTL.chr${chr}_top.pca ${genoCovfile} \
    | awk -v OFS="\t" '\$1=\$1' > Geno10_Pheno10_transcript_Chr${chr}.pca

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
  }' Geno10_Pheno10_transcript_Chr${chr}.pca > xx

  ## Re-order columns according to vcf samples

   awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx sample_listpca.txt  > xy 


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
  }' xy > Geno10_Pheno10_transcript_Chr${chr}.pca

 bgzip Geno10_Pheno10_transcript_Chr${chr}.pca 
  
 awk  'NR==1{for (i=2;i<=NF;i++) print \$i}' pheno_transcript_QTL.chr${chr}.pca > samplefromRNAseq.txt

 
  zcat $genotype | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > vcf_oldVer

 vcftools --vcf vcf_oldVer --keep samplefromRNAseq.txt --maf 0.05 --recode --out Geno_Chr${chr}_transcript

 bgzip Geno_Chr${chr}_transcript.recode.vcf && tabix -p vcf Geno_Chr${chr}_transcript.recode.vcf.gz

 rm xx xy

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
 container 'praveen/qtltools'

 input:
 file (phenotype) 

 tuple val(chr),  file (genotype), file (genoCovfile)

 output:

 
 tuple val(chr), path("pheno_gene_QTL_chr${chr}.bed.*") , emit: phenotype_gene_Bed_ch

 tuple val(chr), file ("Geno10_Pheno10_Geno_Chr${chr}.pca*") , emit:  cis_covariates_gene_ch

 tuple val(chr), path ("Geno_Chr${chr}_gene.recode.vcf.*"),  emit: genotypeQTL_gene_ch 

 
 script:
  //outgeno = genotype[0].toString() - ~/(_UpIds)?(\.vcf)?$/
  //chr = outgeno - ~/ImpuHDktoWGS_RNA_liver_Chr_/

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

 awk -v OFS='\t' 'NR>1 && \$1==${chr} {print}'  Gene_count_Matrices_CommonSamples.tsv | sort -k1,1d -k2,2n -k3,3n > xx

 cat xy xx > pheno_gene_QTL_chr${chr}.bed

 sed -i 's/#Chr/#CHROM/g' pheno_gene_QTL_chr${chr}.bed

 bgzip pheno_gene_QTL_chr${chr}.bed && tabix -p bed pheno_gene_QTL_chr${chr}.bed.gz

 QTLtools pca --bed pheno_gene_QTL_chr${chr}.bed.gz --scale --center --out pheno_gene_QTL_chr${chr}

 awk 'NR<12{ print }' pheno_gene_QTL_chr${chr}.pca > pheno_gene_QTL_chr${chr}_top.pca

 ### Extract common samples and merge common columns from phenotype and genotype covariates

awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} }' pheno_gene_QTL_chr${chr}_top.pca ${genoCovfile} \
   | awk -v OFS="\t" '\$1=\$1' > Geno10_Pheno10_Geno_Chr${chr}.pca

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
  }' Geno10_Pheno10_Geno_Chr${chr}.pca > xx

  ## Re-order columns according to vcf samples

   awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx sample_listpca.txt  > xy 

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
  }' xy > Geno10_Pheno10_Geno_Chr${chr}.pca


 bgzip Geno10_Pheno10_Geno_Chr${chr}.pca 

 awk  'NR==1{for (i=2;i<=NF;i++) print \$i}' pheno_gene_QTL_chr${chr}.pca > samplefromRNAseq.txt

  zcat $genotype | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > vcf_oldVer

 vcftools --vcf vcf_oldVer --keep samplefromRNAseq.txt  --maf 0.05  --recode --out Geno_Chr${chr}_gene

 bgzip Geno_Chr${chr}_gene.recode.vcf && tabix -p vcf Geno_Chr${chr}_gene.recode.vcf.gz

 rm xx xy
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
  container 'praveen/qtltools'
   input:

   file (covariate_s) 

   tuple val(chr), file (phenotype_s),  file (genotype),  file (genoCovfile)

   output:

   tuple val(chr), path ("Geno10_Pheno10_Splice_Chr${chr}.pca.*") , emit: cis_covariates_splice_ch

   tuple val(chr), path ("${phenotype_s}.bed.*") , emit: phenotype_splice_Bed_ch

   tuple val(chr), path  ("Geno_Chr${chr}_splice.recode.*") , emit: genotypeQTL_splice_ch

  script:

  //outgeno = genotype_s[0].toString() - ~/(_UpIds)?(\.vcf)?$/
  //chr = outgeno - ~/ImpuHDktoWGS_RNA_liver_Chr_/

 """

bcftools query -l ${genotype} > samplelist

sed '1iSampleID' samplelist > sample_listpca.txt


 sed 's/comm_//g' $phenotype_s | awk -F'\t' 'NR>1{print \$1,\$2,\$3,\$4}' | awk -F'_' '{print \$1"_"\$2"\t"\$3}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$1":"\$2"-"\$3,\$4,\$5}' | sed '1i#Chr\tstart\tend\tpid\tgid\tstrand' > ${phenotype_s}_MOD.bed.header

 sed 's/comm_//g' $phenotype_s  | awk  '{for(i=5;i<=NF;i++) printf \$i"\t"; print ""}'  > ${phenotype_s}_MOD.bed.counts

 paste -d '\t' ${phenotype_s}_MOD.bed.header ${phenotype_s}_MOD.bed.counts | awk -v OFS='\t' 'NR==1 {print}' > xx

 paste -d'\t' ${phenotype_s}_MOD.bed.header ${phenotype_s}_MOD.bed.counts |  awk -v OFS='\t' 'NR>1 {print}' | sort -k1,1d -k2,2n -k3,3n > xy

 cat xx xy | sed 's/_liver//g' > ${phenotype_s}.bed

 rm ${phenotype_s}_MOD.bed.header ${phenotype_s}_MOD.bed.counts

 bgzip  ${phenotype_s}.bed && tabix -p bed ${phenotype_s}.bed.gz

 sed 's/comm_//g'  $covariate_s | awk 'NR==1;NR>1{\$1="splicePC_"\$1 ;print}' | awk -v OFS=' ' '{\$1=\$1}1' >  ${covariate_s}_MOD.pca

 sed  's/SampleID/id/g' ${genoCovfile} > genCovfile_Mod.txt

  ### Extract common samples and merge common columns from phenotype and genotype covariates

 awk 'FNR==NR && FNR==1{ for(i=1;i<=NF;i++){ b[i]=\$i}; print;  i--; next } \
 FNR==NR{ print; next } FNR!=NR && FNR==1{ for(j=1;j<=NF;j++){ c[\$j]=j};  next } \
  FNR!=NR && FNR>1{ for(k=1;k<=i;k++){ printf("%s%s",\$c[b[k]],k==i?RS:FS)} } '  ${covariate_s}_MOD.pca genCovfile_Mod.txt \
    > Geno10_Pheno10_Splice_Chr${chr}.pca

  awk 'NR==1{for (i=2; i<=NF; i++) print \$i}' ${covariate_s}_MOD.pca > Ind_leafcutter.txt

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
  }' Geno10_Pheno10_Splice_Chr${chr}.pca > xx

  ## Re-order columns according to vcf samples

   awk 'NR==FNR {a[\$1]=\$0;next} \$1 in a {print a[\$1]}' xx sample_listpca.txt  > xy 

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
  }' xy > Geno10_Pheno10_Splice_Chr${chr}.pca



  bgzip Geno10_Pheno10_Splice_Chr${chr}.pca 

  zcat $genotype | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > vcf_oldVer

  vcftools --vcf vcf_oldVer --keep Ind_leafcutter.txt --maf 0.05  --recode --out Geno_Chr${chr}_splice

  bgzip Geno_Chr${chr}_splice.recode.vcf && tabix -p vcf Geno_Chr${chr}_splice.recode.vcf.gz

  rm xx xy

   """
}
