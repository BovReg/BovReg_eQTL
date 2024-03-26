   
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
    file "star" 

    script:
      def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
      """
      mkdir star
      STAR  --runMode genomeGenerate \
            --runThreadN ${task.cpus} \
            --sjdbGTFfile $gtf \
            --genomeDir star/ \
            --genomeFastaFiles $fasta \
            $avail_mem
        """
}