nextflow.enable.dsl=2

/*   ------------ eQTL Nextflow pipeline script02: Extract genotype from samples having corresponding RNAseq samples ----------------- */

/* $projectDir is the current working directory where the script is located */



/* NOTE : Take care of input and output paths before running the script  */


/**  --Input: Input genotype data  **/
/* NOTE: The number of chromsomes can be altered based on user requirements or spcies on interest */
Channel
     .from (25..26)
     .map{chr -> tuple("${chr}",file("$projectDir/Demodata/Demo_genotype/ImpWGS_Ars1.2_Chr${chr}.vcf.gz"))}
     .set {genotype_input_ch}



/**  --Input: Text file with the list of coressponding sample identifiers with genotype and RNAseq data in two columns -- **/ 

params.Corresponding_SampleInfo = "$projectDir/Demodata/RNA_WGS_CorresID.txt"

                                /*  --Output: Path to output files-- */

params.outputGeno = "$projectDir/Output"

                                      def inputfiles() {

                                        log.info """\
                                               Extract genotype data - N F   P I P E L I N E  F O R eQ T L DSL2
                                               ===================================================================
                                               SampleInfo            : ${params.Corresponding_SampleInfo}
                                               Output                : ${params.outputGeno}

                                               """
                                               .stripIndent()
                                       }

                                      inputfiles()


/* Channel objects */
/* single channel objects */

sample_info_ch = Channel.fromPath(params.Corresponding_SampleInfo)




/* calling different process modules */
include {tissuewise_extractGenotype} from './modules_dsl2/Tissuewise_Extract_Merge_Genotypes'



workflow {

 tissuewise_extractGenotype(genotype_input_ch,sample_info_ch.collect())

}

