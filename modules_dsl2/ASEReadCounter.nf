//  Merge the fast fai and dict channels to aSEReadCounter 

process  readCounter{
  tag "$chr"
  container 'broadinstitute/gatk:4.0.8.1'   
  publishDir "${params.outdir}/ASEReadCounter", mode:'copy'


  input:
        file fasta 

        file fastaindex

        file fastadict

        tuple val(chr), file(filtered_vcf), file(fai) 

        file (input_bam) 

    

  output:
     
      tuple val(chr), file ("Chr${chr}_ASE.table.txt")

  script:

      //filename1 = input_bam1[0].toString() - ~/(_SplitNCigarReads)?(\.bam)?$/ 

      //filename2 = input_bam2[0].toString() - ~/(_SplitNCigarReads)?(\.bam)?$/ 

         """
          ## Get the ASE matrix.

          gatk ASEReadCounter \
            -R ${fasta} \
            -I ${input_bam1} \
            -V ${filtered_vcf} \
             -min-depth 10 \
            -O Chr${chr}_ASE.table.txt

         """
}
