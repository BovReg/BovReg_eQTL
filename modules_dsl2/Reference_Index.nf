   
/*
========================================================================================
    STEP 0 - PREPROCESSING - Build STAR index
========================================================================================
*/


process  makeSTARindex {

    tag "star"
    container 'mgibio/star'
    publishDir "${params.starindex}", mode:'copy'

    input:

     file fasta
     file gtf 

    output:
    file "starIndex" 

    script:
      """
      mkdir starIndex
      STAR  --runMode genomeGenerate \
            --runThreadN 10 \
            --sjdbGTFfile $gtf \
            --genomeDir starIndex/ \
            --genomeFastaFiles $fasta \
        """
}