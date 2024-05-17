nextflow.enable.dsl=2

/*
========================================================================================
                                     eQTL-Detect 

   NOTE - Before running the script please take care of the input and output paths                   
========================================================================================
*/

log.info """\
 ==============================================================================================================
                     eQTL-DETECT  [ *** SNP ---> Expression *** ]
 ==============================================================================================================
 starIndex            : ${params.starindex}
 fasta                : ${params.fasta}
 gtf                  : ${params.gtf}
 Input_reads paired   : ${params.pairedreads}
 Input_reads single   : ${params.singlereads}
 SampleInfo           : ${params.corresponding_SampleInfo}
 output               : ${params.outdir}
 CountMatrices        : ${params.countMatrices}
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


Channel.fromPath(params.countMatrices)
    .ifEmpty { error "Cannot find countmatices file in: ${params.countMatrices}" }
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.count_matrices, file(row.gene_count_matrix), file(row.transcript_count_matrix), file(row.splice_count_matrix),file(row.splice_count_pcs))}
    .set { count_matrices_ch }

 // single channel objects
fasta = Channel.fromPath(params.fasta)
gtf = Channel.fromPath(params.gtf)
sample_info_ch = Channel.fromPath(params.corresponding_SampleInfo)


/*
========================================================================================
    Parameters declared in nextflow.config
========================================================================================
*/

anchor_length_ch = params.anchor_length

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum

phenotype_PCs_sQTL_ch = params.phenotype_PCs_sQTL

nominal_cis_ch = params.cis_nominal

permutations_cis_ch = params.cis_permutations

fdr_rate_cis_ch = params.cis_FDR

threshold_trans_ch = params.trans_threshold

permutations_trans_ch = params.trans_permutations

genotype_pcs_ch = params.genotype_pcs

phenotype_PCs_cis = params.phenotype_PCs_cis

phenotype_PCs_trans = params.phenotype_Pcs_trans

/*
========================================================================================
     Modules
========================================================================================
*/

include { makeSTARindex } from './modules_dsl2/Reference_Index'

include {tissuewise_extractGenotype} from './modules_dsl2/Tissuewise_Extract_Merge_Genotypes'

include { genotypeStratificationPCA } from './modules_dsl2/GenotypePCs'



/*
========================================================================================
     Sub-Workflows
========================================================================================
*/


include {PAIREDEND_END_READS} from "./subworkflows/pairedEndReads.nf"

include {SINGLE_END_READS} from "./subworkflows/singleEndReads.nf"

include {QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS} from "./subworkflows/QuantifcationExpressionCounts"

include {ciseQTL_workflow} from "./subworkflows/cis-eQTL_Sub-Wf"

include {transeQTL_workflow} from "./subworkflows/trans-eQTL_Sub-Wf"



workflow {

    main:

    if ( params.countMatrices_input ) {
          
          tissuewise_extractGenotype(genotype_input_ch,sample_info_ch.collect())

          commonSampleIds_ch = tissuewise_extractGenotype.out.genoVcfData_ch
          norm_Ge_count_ch   = count_matrices_ch.map {it[1]}
          norm_Tr_count_ch   = count_matrices_ch.map {it[2]}
          splice_count_ch    = count_matrices_ch.map {it[3]}
          splice_pc_ch       = count_matrices_ch.map {it[4]}

       }

     else {

          makeSTARindex(fasta,gtf)
          
          if(params.pairedEnd_reads){

           PAIREDEND_END_READS(read_pairs_ch, gtf, makeSTARindex.out, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, stringtie_ip_di, leafcutter_ip_di, leafcutter_ip_idx)
           
              stringtie_ip_ch            = PAIREDEND_END_READS.out.stringtie_ip

              leafcutter_bamtojunc_ip_ch = PAIREDEND_END_READS.out.leafcutter_bamtojunc_ip
           }

          if (params.sigleEnd_reads){

           SINGLE_END_READS(single_endReads_ch, gtf, makeSTARindex.out, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, stringtie_ip_di, leafcutter_ip_di, leafcutter_ip_idx)
              
           stringtie_ip_ch            = SINGEL_END_READS.out.stringtie_ip

           leafcutter_bamtojunc_ip_ch = SINGEL_END_READS.out.leafcutter_bamtojunc_ip           

           }

          QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS(stringtie_ip_ch,leafcutter_bamtojunc_ip_ch, gtf, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, genotype_input_ch, sample_info_ch, phenotype_PCs_sQTL_ch)


          commonSampleIds_ch = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.commonSampleIds
          norm_Tr_count_ch   = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.norm_Tr_count
          norm_Ge_count_ch   = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.norm_Ge_count
          splice_pc_ch       = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.splice_pc
          splice_count_ch    = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.splice_Count

     }


    genotypeStratificationPCA(genotype_input_ch, genotype_pcs_ch)

    ciseQTL_workflow(genotypeStratificationPCA.out,commonSampleIds_ch,norm_Tr_count_ch,norm_Ge_count_ch,splice_pc_ch,splice_count_ch,phenotype_PCs_cis,nominal_cis_ch,permutations_cis_ch,fdr_rate_cis_ch)

    transeQTL_workflow(genotypeStratificationPCA.out,commonSampleIds_ch,norm_Tr_count_ch,norm_Ge_count_ch,splice_pc_ch,splice_count_ch,phenotype_PCs_trans,threshold_trans_ch,permutations_trans_ch) 
    
}

                                /* ****END****** */
