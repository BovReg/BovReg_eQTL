 process starAlign1 {
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned_ASE", mode:'copy'

    input:
    tuple val(sample_id), file (bam_pair)
    file index 
    file gtf

    output:
    file "*Log.final.out" 
    file "*.out" 
    file "*SJ.out.tab"
    file "*Log.out" 
    tuple val(sample_id), file ("${sample_id}_ASE_Aligned.sortedByCoord.out.bam"), emit: ase_bam_ch1

    script:

    """

     STAR --genomeDir $index \
             --sjdbGTFfile $gtf \
             --readFilesIn  $bam_pair \
             --twopassMode Basic \
             --alignIntronMin 20 \
             --alignIntronMax  500000 \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${sample_id}_ASE_ \
             --runThreadN 10
        """
}



 process starAlign2 {
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned_ASE", mode:'copy'

    input:
    tuple val(sample_id), file (bam_pair) 
    file index 
    file gtf
    
    output:
    file "*Log.final.out"
    file "*.out" 
    file "*SJ.out.tab"
    file "*Log.out" 
    tuple val(sample_id), file ("${sample_id}_ASE_default_Aligned.sortedByCoord.out.bam"), emit: ase_bam_ch2

    script:

    """

     STAR --genomeDir $index \
             --sjdbGTFfile $gtf \
             --readFilesIn  $bam_pair \
             --twopassMode Basic \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${sample_id}_ASE_default_ \
             --runThreadN 10
        """
}



