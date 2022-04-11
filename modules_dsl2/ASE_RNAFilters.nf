
process addorreplacereadgroup{


  input:
    
      file picardjar 

      tuple val(sample_id), file (align1_bam) 
      tuple val(sample_id), file (align2_bam)

 output:
      tuple val(sample_id), file ("${prefix1}_RG.bam"), emit: bam_rg_ch1
      tuple val(sample_id), file ("${prefix2}_default_RG.bam"), emit: bam_rg_ch2

 script:
      
      prefix1 = align1_bam[0].toString() - ~/(_ASE)?(_Aligned)?(\.sortedByCoord)?(\.out)?(\.bam)?$/
 
      prefix2 = align2_bam[0].toString() - ~/(_ASE)?(_default)?(_Aligned)?(\.sortedByCoord)?(\.out)?(\.bam)?$/

      """
      java -jar $picardjar AddOrReplaceReadGroups \
       I=${align1_bam} \
       O=${prefix1}_RG.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20


       java -jar $picardjar AddOrReplaceReadGroups \
       I=${align2_bam} \
       O=${prefix2}_default_RG.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

      """
}
  




process  markDuplicates {
  
  container 'broadinstitute/gatk:4.0.8.1'   
  publishDir "${params.outdir}/GATK_output_Test1/Picard_MarkDup", mode:'copy'


  input:
      file picardjar 
      tuple val(sample_id), file (rg1_bam) 
      tuple val(sample_id), file (rg2_bam)
        
  output:

      tuple val(sample_id), file ("${prefix1}_dedupped.bam") , emit: markDup_ch1
      tuple val(sample_id), file ("{prefix1}_marked_dup_metrics.txt")
      tuple val(sample_id), file ("${prefix2}_default_dedupped.bam") , emit: markDup_ch2
      tuple val(sample_id), file ("{prefix2}_default_marked_dup_metrics.txt")

  script:

        prefix1 = rg1_bam[0].toString() - ~/(\_RG)?(\.bam)?$/
 
        prefix2 = rg2_bam[0].toString() - ~/(\_default)?(\_RG)?(\.bam)?$/

        """
        java -jar  $picardjar \
        MarkDuplicates I=${rg1_bam} O=${prefix1}_dedupped.bam \
        CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={prefix1}_marked_dup_metrics.txt

        java -jar  $picardjar \
        MarkDuplicates I=${rg2_bam} O=${prefix2}_default_dedupped.bam \
        CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={prefix2}_default_marked_dup_metrics.txt
       """ 
}


//  Merge the fast fai and dict channels to splitNCigarReads_ch /


process splitNCigarReads {
  
  container 'broadinstitute/gatk:4.0.8.1'   
  publishDir "${params.outdir}/GATK_output_Test1/SplitNCigarReads", mode:'copy'


  input:
       file fasta 

       file fastaindex

       file fastadict
       
       tuple val(sample_id), file (picard_op1) 

       tuple val(sample_id), file (picard_op2)

  output:

      tuple val(sample_id), file ("${filename1}_SplitNCigarReads.bam"), emit: cigarReads_ch1

      tuple val(sample_id), file ("${filename2}_SplitNCigarReads.bam"), emit: cigarReads_ch2

  script:

     filename1 = picard_op1[0].toString() - ~/(_dedupped)?(\.bam)?$/

     filename2 = picard_op2[0].toString() - ~/(_dedupped)?(\.bam)?$/

    """
       ##2. Split N trim and reassign mapping qualities
    #Reassign OneMapping Quality: reassign all good alignments (MAPQ=225) to the default value of 60. 
    
    gatk SplitNCigarReads \
        --spark-runner LOCAL \
        -R ${fasta} \
        -I ${picard_op1} \
        -O ${filename1}_SplitNCigarReads.bam

    gatk SplitNCigarReads \
        --spark-runner LOCAL \
        -R ${fasta} \
        -I ${picard_op2} \
        -O ${filename2}_SplitNCigarReads.bam

    """
}

