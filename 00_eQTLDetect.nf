nextflow.enable.dsl=2
/*   ------------------- eQTL Nextflow pipeline script00: Built STAR index -------------------------- */

/* $projectDir is the current working directory */

/**           ****************** NOTE ******************                           */

/*   Take care of input and output paths before running the script                  */



/* Input: Fasta and gtf files of desired species  */

/* $projectDir is the current working directory where the script is located */

params.fasta = "$projectDir/Reference_Genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
params.gtf = "$projectDir/Reference_Genome/Bos_taurus.ARS-UCD1.2.109.gtf"

/* Output: Path to starindex */

params.starindex = "$projectDir/GenomeIndex_Ensembl109"


                                log.info """\
                                         R N A S E Q STAR INDEX BUILT - N F   P I P E L I N E  F O R eQ T L
                                         ==================================================================
                                         fasta        : ${params.fasta}
                                         gtf          : ${params.gtf}
                                         starindex    : ${params.starindex}

                                         """
                                         .stripIndent()

/* Channel objects */

/* single channel objects*/
fasta = file(params.fasta)
gtf = file(params.gtf)

              /* Module calling */
include { makeSTARindex } from './modules_dsl2/Reference_Index'


workflow {

    main:

    makeSTARindex(fasta,gtf)

}
                                /* ****END****** */
