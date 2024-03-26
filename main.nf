nextflow.enable.dsl=2

/*
========================================================================================
                                     eQTL-Detect 

                $projectDir is the current working directory 

   NOTE - Before running the script please take care of the input and output paths                   
========================================================================================
*/


/*
========================================================================================
  Input data (Fasta and gtf files of desired species) 
========================================================================================
*/

params.fasta = "$projectDir/Reference_Genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"

params.gtf = "$projectDir/Reference_Genome/Bos_taurus.ARS-UCD1.2.109.gtf"
/* trimmomatic file to perform trimming */
params.trimmomaticjar = "$projectDir/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"

/* trimming adapter for pair-end reads */
params.adapter_PE = "$projectDir/Softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

/* trimming adapter for single-end reads */
params.adapter_SE = "$projectDir/Softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

/*  1. RNAseq: paired-end library first stranded */
params.reads="$projectDir/Demo_RNAseqFastq_BovReg/*_R{1,2}_subset.fastq.gz"

params.singlereads=""
/*  2. RNAseq: paired-end library unstranded */
//empty

/*  3. RNAseq: paired-end library second stranded */

/**  --Input: Text file with the list of coressponding sample identifiers with genotype and RNAseq data in two columns -- **/ 

params.Corresponding_SampleInfo = "$projectDir/RNA_WGS_CorresID_BovReg.txt"


/* NOTE: The number of chromsomes can be altered based on user requirements or spcies on interest */

starChr_ch = params.starChr

endChr_ch = params.endChr
/* alter the path of genotype file based on user requirement*/
Channel
     .from (starChr_ch..endChr_ch)
     .map{chr -> tuple("${chr}",file("$projectDir/Demo_genotype_BovReg/Bovreg_demogeno_Chr${chr}.vcf.gz"))}
     .set {genotype_input_ch}


/*  --Path to supporting files: QTLtools R script for FDR correction of cis and trans-eQTL hits -- */
params.FDR_cis = "$projectDir/modules_dsl2/runFDR_cis_QTLtools.R"
params.FDR_trans = "$projectDir/modules_dsl2/qtltools_runFDR_ftrans_Mod.R"


/* This will cluster together the introns fond in the junc files listed in junction_files.txt (created below in workflow), 
  requiring 50 split reads supporting each cluster and allowing introns of up to 500kb */
params.leafcutter_cluster = "$projectDir/leafcutter/scripts/leafcutter_cluster_Bovine_regtools.py"

/* This scripts a) calculates intron excision ratios 
  b) filter out introns used in less than 40% of individuals or with almost no variation 
  c) output these ratios as gzipped txt files along with a user-specified number of PCs*/
params.phenotype_table = "$projectDir/leafcutter/scripts/prepare_phenotype_table.py"


/* List of samples common in genotype data. NOTE: The list can be extracted from any chromosome file */

params.sampleIND_gene = "$projectDir/GenoSamples_tsv_Chr2.txt"

params.sampleIND_transcript = "$projectDir/GenoSamples_gtf_Chr2.txt"

params.sampleIND_splice = "$projectDir/GenoSamples_splice_Chr2.txt"


/* Input count matrices and splice PCs*/

params.gene_count_Matrix = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/Count_Matrices_Input/Gene_count_Matrices_filtered.tsv"

params.transcript_count_Matrix = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/Count_Matrices_Input/Transcript_count_Matrices_filtered.tsv"

params.splice_count_Matrix = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/Count_Matrices_Input/Splice_count_Matrices_filtered.tsv"

params.splice_pcs = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/Count_Matrices_Input/Splice_pcs.txt"


/* bam files input data after alignment */

params.stringtie_ip_data = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/samsort_StringTie/*"

params.leafcutter_ip_data = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/samsort_leafcutter/*.bam"

params.leafcutter_ip_idx = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/samsort_leafcutter/*.bai"


/*
========================================================================================
 Output: Path to output files
========================================================================================
*/

params.outdir = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output"

/* Output: Path to starindex */

params.starindex = "$projectDir/../../../disk3/chitneedi/BovReg_eQTL_Output/GenomeIndex_Ensembl109"


                                log.info """\
                                         ===================================================================================
                                                         eQTL-DETECT  (SNP)--->(" *expression* ")
                                         ===================================================================================
                                         fasta        : ${params.fasta}
                                         gtf          : ${params.gtf}
                                         starindex    : ${params.starindex}
                                         Input        : ${params.reads}
                                         gtf          : ${params.gtf}
                                         SampleInfo   : ${params.Corresponding_SampleInfo}
                                         output       : ${params.outdir}
                                         ===================================================================================
                                         """
                                         .stripIndent()
/*
========================================================================================
    Channel objects 
========================================================================================
*/

/* single channel objects*/
fasta = file(params.fasta)
gtf = file(params.gtf)
trimmomaticjar_ch = file(params.trimmomaticjar)
adapter_PE_ch=file(params.adapter_PE)
adapter_SE_ch=file(params.adapter_SE)
ch_FDR_cis = file(params.FDR_cis)
ch_FDR_trans = file(params.FDR_trans)

sample_info_ch = Channel.fromPath(params.Corresponding_SampleInfo)

Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set {gtf_ARS}


     /* Reads with first stranded library */

read_pairs_ch = Channel.fromFilePairs(params.reads)     

/*
Channel
      .fromPath(params.singlereads)
      .ifEmpty { exit 1, " file not found: ${params.singlereads}" }
      .set {single_endReads_ch}
*/

Channel
        .fromPath(params.gene_count_Matrix)
        .ifEmpty { exit 1, " file not found: ${params.gene_count_Matrix}" }
        .set {gene_count_Matrix_ch}

Channel
        .fromPath(params.transcript_count_Matrix)
        .ifEmpty { exit 1, " file not found: ${params.transcript_count_Matrix}" }
        .set {transcript_count_Matrix_ch}

Channel
        .fromPath(params.splice_count_Matrix)
        .ifEmpty { exit 1, " file not found: ${params.splice_count_Matrix}" }
        .set {splice_count_Matrix_ch}

Channel
        .fromPath(params.splice_pcs)
        .ifEmpty { exit 1, " file not found: ${params.splice_pcs}" }
        .set {splice_pcs_ip}

Channel
        .fromPath(params.stringtie_ip_data)
         .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_merged_sorted','').replaceAll('.bam','')
             [ name, content]  }
        .set {stringtie_ip_di}

Channel
        .fromPath(params.leafcutter_ip_data)
         .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_leafcutter_merged_sorted','').replaceAll('.bam','')
             [ name, content]  }
        .set {leafcutter_ip_di}


Channel
        .fromPath(params.leafcutter_ip_idx)
         .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_leafcutter_merged_sorted','').replaceAll('.bai','')
             [ name, content]  }
        .set {leafcutter_ip_idx}

/*
========================================================================================
    parameters declared in json file 
========================================================================================
*/

anchor_length_ch = params.anchor_length

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum

ch_leafcutter_cluster_py = file(params.leafcutter_cluster)

ch_leafcutter_table_py = file(params.phenotype_table)

phenotype_PCs_sQTL_ch = params.phenotype_PCs_sQTL

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum

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
    Include Modules
========================================================================================
*/
include { makeSTARindex } from './modules_dsl2/Reference_Index'

include {tissuewise_extractGenotype} from './modules_dsl2/Tissuewise_Extract_Merge_Genotypes'

include { genotypeStratificationPCA } from './modules_dsl2/GenotypePCs'



/*
========================================================================================
    Include Sub-Workflows
========================================================================================
*/


include {PAIREDEND_END_READS} from "./subworkflows/pairedEndReads.nf"

include {SINGLE_END_READS} from "./subworkflows/singleEndReads.nf"

include {QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS} from "./subworkflows/QuantifcationExpressionCounts"

include {ciseQTL_workflow} from "./subworkflows/cis-eQTL_Sub-Wf"

include {transeQTL_workflow} from "./subworkflows/trans-eQTL_Sub-Wf"

params.countMatrices_input = false

params.pairedEnd_reads = false

params.sigleEnd_reads = false



workflow {

    main:

    if ( params.countMatrices_input ) {

          tissuewise_extractGenotype(genotype_input_ch,sample_info_ch.collect())

          commonSampleIds_ch = tissuewise_extractGenotype.out.genoVcfData_ch
          norm_Tr_count_ch   = transcript_count_Matrix_ch
          norm_Ge_count_ch   = gene_count_Matrix_ch
          splice_pc_ch       = splice_pcs_ip
          splice_count_ch    = splice_count_Matrix_ch

       }

     else {

          makeSTARindex(fasta,gtf)
          
          if(params.pairedEnd_reads){

           PAIREDEND_END_READS(read_pairs_ch, trimmomaticjar_ch, adapter_PE_ch, gtf_ARS, makeSTARindex.out, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, stringtie_ip_di, leafcutter_ip_di, leafcutter_ip_idx)
           
              stringtie_ip_ch            = PAIREDEND_END_READS.out.stringtie_ip

              leafcutter_bamtojunc_ip_ch = PAIREDEND_END_READS.out.leafcutter_bamtojunc_ip
           }

          if (params.sigleEnd_reads){

           SINGEL_END_READS(single_endReads_ch, trimmomaticjar_ch, adapter_SE_ch, gtf_ARS, makeSTARindex.out, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, stringtie_ip_di, leafcutter_ip_di, leafcutter_ip_idx)
              
           stringtie_ip_ch            = SINGEL_END_READS.out.stringtie_ip

           leafcutter_bamtojunc_ip_ch = SINGEL_END_READS.out.leafcutter_bamtojunc_ip           

           }

          QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS(stringtie_ip_ch,leafcutter_bamtojunc_ip_ch, gtf_ARS, anchor_length_ch, intron_length_minimum_ch, intron_length_maximum_ch, genotype_input_ch, sample_info_ch, phenotype_PCs_sQTL_ch, ch_leafcutter_cluster_py, ch_leafcutter_table_py)


          commonSampleIds_ch = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.commonSampleIds
          norm_Tr_count_ch   = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.norm_Tr_count
          norm_Ge_count_ch   = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.norm_Ge_count
          splice_pc_ch       = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.splice_pc
          splice_count_ch    = QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS.out.splice_Count

     }

    genotypeStratificationPCA(genotype_input_ch, genotype_pcs_ch)

    ciseQTL_workflow(genotypeStratificationPCA.out,commonSampleIds_ch,norm_Tr_count_ch,norm_Ge_count_ch,splice_pc_ch,splice_count_ch,ch_FDR_cis,phenotype_PCs_cis,nominal_cis_ch,permutations_cis_ch,fdr_rate_cis_ch)

    transeQTL_workflow(genotypeStratificationPCA.out,commonSampleIds_ch,norm_Tr_count_ch,norm_Ge_count_ch,splice_pc_ch,splice_count_ch,ch_FDR_trans,phenotype_PCs_trans,threshold_trans_ch,permutations_trans_ch) 
    
}

                                /* ****END****** */
