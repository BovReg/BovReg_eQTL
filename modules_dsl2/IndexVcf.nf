
process indexvcf_ASE{
     tag "$chr"
     container 'broadinstitute/gatk:4.0.8.1'   
     publishDir "${params.outdir}/ASEReadCounter", mode:'copy'

     input:

     tuple val(chr), file(genotype_input_ch)

     output:
     
     tuple  val(chr),file ("${genotype_input_ch}"), emit: vcf 
     tuple  val(chr), file ("${genotype_input_ch}.tbi"), emit: vcf_index 
     script:
      
        //prefix = vcf[0].toString() 
        //name = prefix - ~/ImpuHDktoWGS_segfam_ChromMas_NUDA_Ars1.2_NdDidierfile_1043IND_/ 

        """
            ## Index vcf file

              gatk IndexFeatureFile \
                -F ${genotype_input_ch} 

              mv ${genotype_input_ch} temp

              mv temp ${genotype_input_ch}

           """

}
