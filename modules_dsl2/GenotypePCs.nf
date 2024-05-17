
/* 

process GenotypeStratificationPCA: Perform pca analysis using R package SNPRelate for genotype data.
   - The input is the imputed WGS performed per chromosome and top 10 PCs were filtered.
   - The top 10 PCs were forwared to the processes : TranscriptCountsPCA,RNAsplicePCS and GeneCountsPCA.
   - In these processes these genotype PCs were combined with PCs from RNAseq at transcript, splce and gene level per chromsome.

*/

process genotypeStratificationPCA {
 tag "on chromosome ${chr}"
 publishDir "${params.outdir}/GenotypeCov", mode:'copy'
 container 'praveenchitneedi/rpackages-eqtl:1.0.1' 


  input:
    tuple val(chr), file(vcf_input)

    val(genotype_pcs)

  output:

    tuple val(chr), file ("GenoCovar${chr}.txt") 

  script:

     """
     #!/usr/bin/env Rscript

     library(SNPRelate)
     library(dplyr)

     vcf.fn<-paste0("$vcf_input")

     snpgdsVCF2GDS(vcf.fn,"ccm.gds",method="biallelic.only")
     genofile<-snpgdsOpen("ccm.gds")
     ccm_pca<-snpgdsPCA(genofile)
     ### Plot to check for the variantion explained by each PCA
       # plot(ccm_pca,eig=c(1L,2L,3L,4L))
     ### Extract eigen vectors 
     pca_genotype<-ccm_pca\$eigenvect[,1:$genotype_pcs]
     rownames(pca_genotype) <- ccm_pca\$sample.id
     pcat <- t(pca_genotype)
     ## convert double to dataframe object
     pcatdf<- as.data.frame(pcat)
     ## Round the decimal points to 5 
     fpcatdf <- pcatdf %>% mutate_if(is.numeric, round, digits = 5)
     write.table(data.frame("SampleID"=rownames(fpcatdf),fpcatdf),"GenoCovar${chr}.txt", quote=FALSE,row.names=FALSE)
     """
}