//  Merge the fast fai and dict channels to aSEReadCounter 

process  readCounter{
  
  container 'broadinstitute/gatk:4.0.8.1'   
  publishDir "${params.outdir}/ASEReadCounter_STARWASP", mode:'copy'


  input:
        file fasta 

        file fastaindex

        file fastadict

        tuple val(name), file(filtered_vcf), file(fai) 

        tuple val(sample_id), file (input_bam1) 

        tuple val(sample_id), file (input_bam2) 

  output:
     
       tuple val(sample_id),file ("${filename1}.ASE.table_Cigar.txt")

       tuple val(sample_id),file ("${filename2}.ASE.table_Cigar.txt")

  script:

      filename1 = input_bam1[0].toString() - ~/(_SplitNCigarReads)?(\.bam)?$/ 

      filename2 = input_bam2[0].toString() - ~/(_SplitNCigarReads)?(\.bam)?$/ 

         """
           ## Get the ASE matrix.

          gatk ASEReadCounter \
            -R ${fasta} \
            -I ${input_bam1} \
            -V ${filtered_vcf} \
             -min-depth 10 \
            -O ${filename1}.ASE.table_Cigar.txt


          gatk ASEReadCounter \
            -R ${fasta} \
            -I ${input_bam2} \
            -V ${filtered_vcf} \
            -min-depth 10 \
            -O ${filename2}.ASE.table_Cigar.txt  

         """
}
