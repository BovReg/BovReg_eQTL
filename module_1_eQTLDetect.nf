nextflow.enable.dsl=2

/*  

========================================================================================
                                   eQTL-Detect

             Module-01: RNAseq read trimming, alignment and quantification

 NOTE : Take care of input and output paths before running the script  

========================================================================================
*/



log.info """\

==============================================================================================================
                     eQTL-DETECT  [ *** SNP ---> Expression *** ] script 01
 ==============================================================================================================
 starIndex            : ${params.starindex}
 fasta                : ${params.fasta}
 gtf                  : ${params.gtf}
 Input_reads paired   : ${params.pairedreads}
 Input_reads single   : ${params.singlereads}
 SampleInfo           : ${params.corresponding_SampleInfo}
 output               : ${params.outdir}
 Bam_input_files      : ${params.bamIpfiles}
 Geno_input_files     : ${params.genoIpfiles}
 ==============================================================================================================

"""
.stripIndent()



/*
========================================================================================
    Channel objects declared in nextflow.config
========================================================================================
*/
    // Channels from tsv files
Channel.fromPath(params.genoIpfiles)
        .ifEmpty { error "Cannot find countmatices file in: ${params.genoIpfiles}" }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.chromosome, file(row.vcfFile))  }
        .set {genotype_input_ch}

Channel.fromPath(params.pairedreads)
    .ifEmpty { error "Cannot find countmatices file in: ${params.pairedreads}" }
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.sampleId, file(row.read1), file(row.read2))}
    .set { read_pairs_ch }

Channel.fromPath(params.bamIpfiles)
        .ifEmpty { error "Cannot find countmatices file in: ${params.bamIpfiles}" }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sampleId, file(row.stringTieBam))  }
        .set {stringtie_ip_di}


Channel.fromPath(params.bamIpfiles)
        .ifEmpty { error "Cannot find countmatices file in: ${params.bamIpfiles}" }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sampleId, file(row.leafcutterBam))  }
        .set {leafcutter_ip_di}


Channel.fromPath(params.bamIpfiles)
        .ifEmpty { error "Cannot find countmatices file in: ${params.bamIpfiles}" }
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sampleId, file(row.leafcutterBai))  }
        .set {leafcutter_ip_idx}

/* single channel objects*/
fasta = Channel.fromPath(params.fasta)
gtf = Channel.fromPath(params.gtf)

  
/*
========================================================================================
    Parameters declared in nextflow.config
========================================================================================
*/

anchor_length_ch = params.anchor_length

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum

phenotype_PCs_sQTL_ch = params.phenotype_PCs_sQTL

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum


/*
========================================================================================
     Modules
========================================================================================
*/

include { makeSTARindex } from './modules_dsl2/Reference_Index'


/*
========================================================================================
     Sub-Workflows
========================================================================================
*/


include {PAIREDEND_END_READS} from "./subworkflows/pairedEndReads.nf"

include {SINGLE_END_READS} from "./subworkflows/singleEndReads.nf"
                             

/*
========================================================================================
     Sub-Workflows
========================================================================================
*/


include {PAIREDEND_END_READS} from "./subworkflows/pairedEndReads.nf"

include {SINGLE_END_READS} from "./subworkflows/singleEndReads.nf"


workflow {

    main:

        makeSTARindex(fasta,gtf)

          if(params.pairedEnd_reads){

           PAIREDEND_END_READS(read_pairs_ch, gtf, makeSTARindex.out, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, stringtie_ip_di, leafcutter_ip_di, leafcutter_ip_idx)
           
           }

          if (params.sigleEnd_reads){

           SINGLE_END_READS(single_endReads_ch, gtf, makeSTARindex.out, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, stringtie_ip_di, leafcutter_ip_di, leafcutter_ip_idx)
                        

           } 
    
}





                                              /* ****END****** */




