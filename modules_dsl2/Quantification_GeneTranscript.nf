                           /*
                             * STEP 5 - Stringtie (FPKM quantification, reads and transcripts)
                             */


process stringtieQuant_library_firststranded {
    tag "stringtieQuant on $sample_id"
    container 'gawbul/docker-stringtie'
    publishDir "${params.outdir}/stringtieQuantify", mode:'copy'
    input:

    tuple val(sample_id), file (star_sorted)

    file(gtf) 

    output:

    //tuple  file ("*.gtf") , file ("*.tsv") 
    
    tuple  val(sample_id), file ("*.gtf") , file ("*.tsv")

    script:
  
   """
   stringtie --rf  -p 6 -e -B -G $gtf \
    -o ${sample_id}_stringtie.gtf \
    -A ${sample_id}_stringtie.tsv \
    $star_sorted

  """
}






 //stringtie_bam_ip = stringtie_bam_ch.merge(stringtie_bai_ch)
process stringtieQuant_library_unstranded {
    tag "stringtieQuant on $sample_id"
    container 'gawbul/docker-stringtie'
    publishDir "${params.outdir}/stringtieQuantify", mode:'copy'
    input:

    tuple val(sample_id), file (star_sorted_bam)

    file(gtf) 

    output:
    
    tuple val(sample_id), file ("*.gtf") , file ("*.tsv")

    script:


   """
   stringtie  -p 6 -e -B -G $gtf \
       -o ${sample_id}_stringtie.gtf \
       -A ${sample_id}_stringtie.tsv \
       $star_sorted_bam
  """
}



 //stringtie_bam_ip = stringtie_bam_ch.merge(stringtie_bai_ch)
process stringtieQuant_library_secondstranded {
    tag "stringtieQuant on $sample_id"
    container 'gawbul/docker-stringtie'
    publishDir "${params.outdir}/stringtieQuantify", mode:'copy'
    input:
    //file bam from stringtie_ch
    //tuple val(sample_id), file (star_sorted_bam)

    tuple val(sample_id), file (star_sorted)

    file(gtf) 

    output:

    //tuple  file ("*.gtf") , file ("*.tsv") 
    
    tuple  val(sample_id), file ("*.gtf") , file ("*.tsv")

    script:

  //samplename = star_sorted[0].toString() - ~/(\.bam)?$/
   //name = star_sorted[0].toString()
   //prefix = name.split('_')[0]
  
   """
   stringtie --fr  -p 6 -e -B -G $gtf \
    -o ${sample_id}_stringtie.gtf \
    -A ${sample_id}_stringtie.tsv \
    $star_sorted

  """
}
