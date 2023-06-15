
                                /*
                                 * STEP 1 - Check read quality FastQC
                                 */

process fastqc_PE {
    tag "FASTQC on $prefix"
    container 'nextflow/rnaseq-nf'
    publishDir "${params.outdir}/fastqc_results", mode:'copy'
    input:
    tuple val(prefix), path(pairedEndread_files) 

    output:
     tuple val(prefix), file("fastqc_${prefix}_logs") 

    script:
     name = pairedEndread_files[0].toString()
     prefix = name.split('_')[0]

    """
   
    mkdir fastqc_${prefix}_logs
    fastqc -o fastqc_${prefix}_logs -f fastq -q ${pairedEndread_files}
    """
}



process fastqc_SE {
    tag "FASTQC on $prefix"
    container 'nextflow/rnaseq-nf'
    publishDir "${params.outdir}/fastqc_results", mode:'copy'
    input:
      file(singelEndread_files) 

    output:
     tuple val(prefix), file("fastqc_${prefix}_logs") 


    script:

   name = singelEndread_files[0].toString()
   prefix = name.split('_')[0]


    """
    mkdir fastqc_${prefix}_logs
    fastqc -o fastqc_${prefix}_logs -f fastq -q ${singelEndread_files}
    """
}
                                    /*
                                     * STEP 2 - Read trimming
                                     */


process trimmomatic{
  tag "trimmomatic on $prefix"
  publishDir "${params.outdir}/Trimmomatic", mode:'copy'

  input:
      tuple val(prefix), path(readPair_files) 
      file trimmjar  
      file adapter 
  output:
     tuple val(prefix), file ("${prefix}_R1.clean.fastq.gz"), emit: forward_strand_trim
     tuple val(prefix), file ("${prefix}_R2.clean.fastq.gz"), emit: revese_strand_trim
     tuple val(prefix), file ("${prefix}_R1_unpaired.fastq.gz"), emit: forward_unpaired_trim
     tuple val(prefix), file ("${prefix}_R2_unpaired.fastq.gz"), emit: revese_unpaired_trim
     tuple val(prefix), file ("${prefix}_trimomatic.log")

  script:
   //prefix = reads[0].toString() - ~/(_R1)?(\.fastq)?(\.gz)?$/
   //sample = prefix.split('_')[0]

    name = readPair_files[0].toString()
    prefix = name.split('_')[0]


  """
   java -jar ${trimmjar} PE -phred33 ${readPair_files} ${prefix}_R1.clean.fastq.gz ${prefix}_R1_unpaired.fastq.gz ${prefix}_R2.clean.fastq.gz ${prefix}_R2_unpaired.fastq.gz -threads 10 \
   ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> ${prefix}_trimomatic.log
  """
}



// Trimmomatic report extracted using the following command
//grep -e "Input Read Pairs:" $baseDir/results_eQTL_MappingQunatif/Trimmoatic_output/*_trimomatic.log | awk '{print $1,$4,$7$8,$12$13,$17$18,$20$21}' | sed 's/_trimomatic.log:Input//g' | sed '1iSample Input_Read_Pairs Both_Surviving Forward_Only_Surviving Reverse_Only_Surviving Dropped' | tr ' ' '\t'  >> Trimmomatic_output_report.txt



process trimmomatic_singleEnd{
  tag "trimmomatic on $prefix"
  publishDir "${params.outdir}/Trimmomatic_SE", mode:'copy'

  input:
      file(singelEndread_files) 
      file trimmjar  
      file adapter 
  output:
     tuple val(prefix), file ("${prefix}_SingleEnd.clean.fastq.gz"), emit: forward_strand_trim
     tuple val(prefix), file ("${prefix}_trimomatic.log")

  script:   
   name = singelEndread_files[0].toString()
   prefix = name.split('_')[0]

  """
   java -jar ${trimmjar} SE -phred33 ${singelEndread_files} ${prefix}_SingleEnd.clean.fastq.gz -threads 10 \
   ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> ${prefix}_trimomatic.log
  """
}
