
   /*
                                  STEP 3 - Read alignment with STAR and the output BAM files are supplied
                                            to Stringtie , features counts and leafcutter for read quantification
                                            */


 process starAlign_GeneTranscript{
    tag "starAlingment on $sample_id"
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned", mode:'copy'

    input:
    tuple val(sample_id), file(read1),file(read2)
    file gtf 
    file index
    //set sample_id, file (read_splicevar) from read_splicevar_ch
    output:
    file "*Log.final.out" 
    tuple val(sample_id), file ("${sample_id}_paired_Aligned.sortedByCoord.out.bam"), emit:tomerge_paired_ch
    //file "${sample_id}_Aligned.toTranscriptome.out.bam" into star_Transc_bam
    file "*.out" 
    file "*SJ.out.tab"
    file "*Log.out" 

    script:


    """
     STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn ${read1} ${read2} \
             --runThreadN 10 \
             --outSAMtype BAM SortedByCoordinate \
             --readFilesCommand zcat \
             --outFileNamePrefix ${sample_id}_paired_

    """
}


 process starAlign_Splicing {
    tag "starAlingment on $sample_id"
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned", mode:'copy'

    input:
    tuple val(sample_id), file(read1),file(read2)
    file gtf 
    file index
    //set sample_id, file (read_splicevar) from read_splicevar_ch
    output:
    file "*Log.final.out" 
    file "*.out" 
    file "*SJ.out.tab"
    file "*Log.out" 
    tuple val(sample_id),file ("${sample_id}_paired_leafcutter_Aligned.out.bam"), emit:tomerge_paired_Leaf_ch

    script:
    //prefix = read_splicevar[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

     //prefix = read1[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

    // STEP 3.1  First STAR alignment of RNA and bam file preperation with sorted for stringtie and featurecounts quantification 
    // STEP 3.2 Second STAR alignment of RNA and bam file preperation with unsorted for leafcutter quantification
    // The unpaired reads from trimmomatic are also aligned 

    """
      
     STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn  ${read1} ${read2} \
             --twopassMode Basic \
             --outSAMstrandField intronMotif \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix ${sample_id}_paired_leafcutter_ \
             --runThreadN 10

    """
}


 process starAlign_unpaired {
    tag "starAlingment on $sample_id"
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned", mode:'copy'

    input:
    tuple val(sample_id),file(unpaired_read1),file(unpaired_read2)
    file gtf 
    file index
    //set sample_id, file (read_splicevar) from read_splicevar_ch
    output:
    //file "*Log.final.out" 
    tuple val(sample_id), file ("${sample_id}_unpaired_R1_Aligned.sortedByCoord.out.bam"), emit: tomerge_unpairedR1_ch
    tuple val(sample_id), file ("${sample_id}_unpaired_R2_Aligned.sortedByCoord.out.bam"), emit: tomerge_unpairedR2_ch
 
    tuple val(sample_id),file ("${sample_id}_leafcutter_unpaired_R1_Aligned.out.bam"), emit: tomerge_unpairedLeafR1_ch
    tuple val(sample_id),file ("${sample_id}_leafcutter_unpaired_R2_Aligned.out.bam"), emit: tomerge_unpairedLeafR2_ch

    script:
    //prefix = read_splicevar[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

     //prefix = read1[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

    // STEP 3.1  First STAR alignment of RNA and bam file preperation with sorted for stringtie and featurecounts quantification 
    // STEP 3.2 Second STAR alignment of RNA and bam file preperation with unsorted for leafcutter quantification
    // The unpaired reads from trimmomatic are also aligned 

    """
      STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn ${unpaired_read1}  \
             --runThreadN 10 \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesCommand zcat \
             --outFileNamePrefix ${sample_id}_unpaired_R1_


       STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
            --readFilesIn ${unpaired_read2}  \
             --runThreadN 10 \
             --outSAMtype BAM SortedByCoordinate \
            --readFilesCommand zcat \
            --outFileNamePrefix ${sample_id}_unpaired_R2_
  
        STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn  ${unpaired_read1} \
             --twopassMode Basic \
             --outSAMstrandField intronMotif \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix ${sample_id}_leafcutter_unpaired_R1_ \
             --runThreadN 10
        

        STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn  ${unpaired_read2} \
             --twopassMode Basic \
             --outSAMstrandField intronMotif \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix ${sample_id}_leafcutter_unpaired_R2_ \
             --runThreadN 10
    """
}

 process starAlign_singleEnd_GeneTranscript{
    tag "starAlingment on $sample_id"
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned", mode:'copy'

    input:
    tuple val(sample_id), file(read1)
    file gtf 
    file index
    //set sample_id, file (read_splicevar) from read_splicevar_ch
    output:
    file "*Log.final.out" 
    tuple val(sample_id), file ("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: stringtie_ch
    //file "${sample_id}_Aligned.toTranscriptome.out.bam" into star_Transc_bam
    file "*.out" 
    file "*SJ.out.tab"
    file "*Log.out" 

    script:
    //prefix = read_splicevar[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

     //prefix = read1[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

    // STEP 3.1  First STAR alignment of RNA and bam file preperation with sorted for stringtie and featurecounts quantification 
    // STEP 3.2 Second STAR alignment of RNA and bam file preperation with unsorted for leafcutter quantification

    """
     STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn ${read1}  \
             --runThreadN 10 \
             --outSAMtype BAM SortedByCoordinate \
             --readFilesCommand zcat \
             --outFileNamePrefix ${sample_id}_
    """
}


process starAlign_singleEnd_Splicing{
    tag "starAlingment on $sample_id"
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned", mode:'copy'

    input:
    tuple val(sample_id), file(read1)
    file gtf 
    file index
    //set sample_id, file (read_splicevar) from read_splicevar_ch
    output:
    file "*Log.final.out" 
    //file "${sample_id}_Aligned.toTranscriptome.out.bam" into star_Transc_bam
    file "*.out" 
    file "*SJ.out.tab"
    file "*Log.out" 
    tuple val(sample_id),file ("${sample_id}_leafcutter_Aligned.out.bam"), emit: samindex_ch

    script:
    //prefix = read_splicevar[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

     //prefix = read1[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

    // STEP 3.1  First STAR alignment of RNA and bam file preperation with sorted for stringtie and featurecounts quantification 
    // STEP 3.2 Second STAR alignment of RNA and bam file preperation with unsorted for leafcutter quantification

    """
      
     STAR --genomeDir ${index} \
             --sjdbGTFfile ${gtf} \
             --readFilesIn  ${read1}  \
             --twopassMode Basic \
             --outSAMstrandField intronMotif \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix ${sample_id}_leafcutter_ \
             --runThreadN 10
    """
}