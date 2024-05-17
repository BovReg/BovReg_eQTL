nextflow.enable.dsl=2

/*
========================================================================================
                                     eQTL-Detect
            script03: Performs cis, trans and sQTL mapping   

   NOTE : Before running the script please take care of the input and output paths                   
========================================================================================
*/ 


log.info """\
 ==============================================================================================================
                     eQTL-DETECT  [ *** SNP ---> Expression *** ]
 ==============================================================================================================
 output               : ${params.outdir}
 SampleInfo           : ${params.corresponding_SampleInfo}
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


Channel.fromPath(params.countMatrices)
    .ifEmpty { error "Cannot find countmatices file in: ${params.countMatrices}" }
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.count_matrices, file(row.gene_count_matrix), file(row.transcript_count_matrix), file(row.splice_count_matrix),file(row.splice_count_pcs))}
    .set { count_matrices_ch }

/* single channel objects */

sample_info_ch = Channel.fromPath(params.corresponding_SampleInfo)



/*
========================================================================================
    Parameters declared in nextflow.config
========================================================================================
*/

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

include {tissuewise_extractGenotype} from './modules_dsl2/Tissuewise_Extract_Merge_Genotypes'

include { genotypeStratificationPCA } from './modules_dsl2/GenotypePCs'


/*
========================================================================================
     Sub-Workflows
========================================================================================
*/

include {ciseQTL_workflow} from "./subworkflows/cis-eQTL_Sub-Wf"

include {transeQTL_workflow} from "./subworkflows/trans-eQTL_Sub-Wf"



workflow {

    main:
       
          tissuewise_extractGenotype(genotype_input_ch,sample_info_ch.collect())

          commonSampleIds_ch = tissuewise_extractGenotype.out.genoVcfData_ch
          norm_Ge_count_ch   = count_matrices_ch.map {it[1]}
          norm_Tr_count_ch   = count_matrices_ch.map {it[2]}
          splice_count_ch    = count_matrices_ch.map {it[3]}
          splice_pc_ch       = count_matrices_ch.map {it[4]}

          genotypeStratificationPCA(genotype_input_ch, genotype_pcs_ch)

          ciseQTL_workflow(genotypeStratificationPCA.out,commonSampleIds_ch,norm_Tr_count_ch,norm_Ge_count_ch,splice_pc_ch,splice_count_ch,phenotype_PCs_cis,nominal_cis_ch,permutations_cis_ch,fdr_rate_cis_ch)

          transeQTL_workflow(genotypeStratificationPCA.out,commonSampleIds_ch,norm_Tr_count_ch,norm_Ge_count_ch,splice_pc_ch,splice_count_ch,phenotype_PCs_trans,threshold_trans_ch,permutations_trans_ch)

}