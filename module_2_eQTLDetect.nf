nextflow.enable.dsl=2

/*
========================================================================================
                                   eQTL-Detect

    Module-02: 1. Extract genotype from samples having corresponding RNAseq samples 
               2. Quantification and merging the read and transcript counts generated from stringtie and
                cluster introns found in junction files and estimate covriates for splicing 
                sites based on PCs

 NOTE : Take care of input and output paths before running the script  

========================================================================================
*/



log.info """\
 ==============================================================================================================
                     eQTL-DETECT  [ *** SNP ---> Expression *** ]
 ==============================================================================================================
 SampleInfo           : ${params.corresponding_SampleInfo}
 output               : ${params.outdir}
 CountMatrices        : ${params.countMatrices}
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


// single channel objects

sample_info_ch = Channel.fromPath(params.corresponding_SampleInfo)

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



/*
========================================================================================
     Modules
========================================================================================
*/

/* calling different process modules */
include {tissuewise_extractGenotype} from './modules_dsl2/Tissuewise_Extract_Merge_Genotypes'

/*
========================================================================================
     Sub-Workflows
========================================================================================
*/


include {QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS} from "./subworkflows/QuantifcationExpressionCounts"



workflow {

    main:
 
    tissuewise_extractGenotype(genotype_input_ch,sample_info_ch.collect())

    leafcutter_bamtojunc_ip = leafcutter_ip_di.join(leafcutter_ip_idx)

    QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS(stringtie_ip_di , leafcutter_bamtojunc_ip , gtf, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, genotype_input_ch, sample_info_ch, phenotype_PCs_sQTL_ch)

}

