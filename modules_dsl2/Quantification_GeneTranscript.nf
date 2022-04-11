                           /*
                             * STEP 5 - Stringtie (FPKM quantification, reads and transcripts)
                             */




 //stringtie_bam_ip = stringtie_bam_ch.merge(stringtie_bai_ch)
process stringtieQuant_with_rf {
    //tag "stringtieQuant on $sample_id"
    container 'gawbul/docker-stringtie'
    publishDir "${params.outdir}/stringtieQuantify", mode:'copy'
    input:
    //file bam from stringtie_ch
    //tuple val(sample_id), file (star_sorted_bam)

    tuple val(sample_id), file (star_sorted)

    file(gtf) 

    output:

    //tuple  file ("*.gtf") , file ("*.tsv") 
    
    tuple  file ("*.gtf") , file ("*.tsv")

    script:

  //samplename = star_sorted[0].toString() - ~/(\.bam)?$/
   //name = star_sorted[0].toString()
   //prefix = name.split('_')[0]
  
   """
   stringtie --rf  -p 6 -e -B -G $gtf \
    -o ${sample_id}_stringtie.gtf \
    -A ${sample_id}_stringtie.tsv \
    $star_sorted

  """
}

 //stringtie_bam_ip = stringtie_bam_ch.merge(stringtie_bai_ch)
process stringtieQuant_without_rf {
    //tag "stringtieQuant on $sample_id"
    container 'gawbul/docker-stringtie'
    publishDir "${params.outdir}/stringtieQuantify", mode:'copy'
    input:
    //file bam from stringtie_ch
    //tuple val(sample_id), file (star_sorted_bam)

    tuple val(sample_id), file (star_sorted_bam)

    file(gtf) 

    output:

    //tuple  file ("*.gtf") , file ("*.tsv") 
    
    tuple file ("*.gtf") , file ("*.tsv")

    script:

  //samplename = star_sorted_bam[0].toString() - ~/(_merged)?(_sorted)?(\.bam)?$/

   """
   stringtie  -p 6 -e -B -G $gtf \
       -o ${sample_id}_stringtie.gtf \
       -A ${sample_id}_stringtie.tsv \
       $star_sorted_bam
  """
}

 //stringtie_bam_ip = stringtie_bam_ch.merge(stringtie_bai_ch)
process stringtieQuant {
    //tag "stringtieQuant on $sample_id"
    container 'gawbul/docker-stringtie'
    publishDir "${params.outdir}/stringtieQuantify", mode:'copy'
    input:
    //file bam from stringtie_ch
    //tuple val(sample_id), file (star_sorted_bam)

    tuple val(sample_id), file (star_sorted)

    file(gtf) 

    output:

    //tuple  file ("*.gtf") , file ("*.tsv") 
    
    tuple  file ("*.gtf") , file ("*.tsv")

    script:

  //samplename = star_sorted[0].toString() - ~/(\.bam)?$/
   //name = star_sorted[0].toString()
   //prefix = name.split('_')[0]
  
   """
   stringtie --rf  -p 6 -e -B -G $gtf \
    -o ${sample_id}_stringtie.gtf \
    -A ${sample_id}_stringtie.tsv \
    $star_sorted

  """
}