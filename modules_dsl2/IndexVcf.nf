
process indexvcf_ASE{

     container 'broadinstitute/gatk:4.0.8.1'   
     publishDir "${params.outdir}/ASEReadCounter_STARWASP", mode:'copy'

     input:

    file vcf 

     output:
     
     tuple  val(name),file ("${vcf}"), emit: vcf 
     tuple  val(name), file ("${vcf}.tbi"), emit: vcf_index 
     script:
      
        prefix = vcf[0].toString() 
        name = prefix - ~/ImpuHDktoWGS_segfam_ChromMas_NUDA_Ars1.2_NdDidierfile_1043IND_/ 

        """
            ## Index vcf file

              gatk IndexFeatureFile \
                -F ${vcf} 

              mv ${vcf} temp

              mv temp ${vcf}

           """

}
