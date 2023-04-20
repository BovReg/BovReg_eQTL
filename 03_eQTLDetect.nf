nextflow.enable.dsl=2

/*   ------------ eQTL Nextflow pipeline script03: merging the read and transcript counts generated from stringtie and cluster introns found in junction files estimate covriates for splicing sites based on PCs----------------- */


/**           ****************** NOTE ******************                           */
/* Take care of input and output paths before running the script                  */



 /* $projectDir is the current working directory /home/chitneedi/Disk2chitneedi/NextFlow */
/*  Script parameters */



/**  --Input: Across sample gene and transcript level expression counts and junction files for splicing QTL-- **/  

params.countsGene="$projectDir/Output/stringtieQuantify/*.tsv"

params.countsTranscripts="$projectDir/Output/stringtieQuantify/*.gtf"

params.junc = "$projectDir/Output/leafcutter_bamtojunc/*.junc"


/* List of samples common in genotype data. NOTE: The list can be extracted from any chromosome file */

params.sampleIND_gene = "$projectDir/Output/eQTLGenoData/GenoSamples_tsv_Chr26.txt"

params.sampleIND_transcript = "$projectDir/Output/eQTLGenoData/GenoSamples_gtf_Chr26.txt"

params.sampleIND_splice="$projectDir/Output/eQTLGenoData/GenoSamples_splice_Chr26.txt"


                                  /*  --Output: Path to output files-- */

params.outdir = "$projectDir/Output"

   
                           /*  --Path to supporting files: python scripts to generate sQTL input -  */

/* This will cluster together the introns fond in the junc files listed in junction_files.txt (created below in workflow), 
  requiring 50 split reads supporting each cluster and allowing introns of up to 500kb */
params.leafcutter_cluster = "$projectDir/leafcutter/scripts/leafcutter_cluster_Bovine_regtools.py"

/* This scripts a) calculates intron excision ratios 
  b) filter out introns used in less than 40% of individuals or with almost no variation 
  c) output these ratios as gzipped txt files along with a user-specified number of PCs*/
params.phenotype_table = "$projectDir/leafcutter/scripts/prepare_phenotype_table.py"


                                  def inputfiles() {

                                    log.info """\
                                           Count matrices - N F   P I P E L I N E  F O R eQ T L DSL2
                                           ===================================================================
                                           GeneMatrix            : ${params.countsGene}
                                           TranscriptMatrix      : ${params.countsTranscripts}
                                           OutputDirectory       : ${params.outdir}
                                           """
                                           .stripIndent()
                                   }

                                  inputfiles()



/* Channel objects */

/* single channel objects*/


junc_liver_ch=Channel.fromPath(params.junc)

gene_St_counts_ch = Channel.fromPath(params.countsGene)

transcript_counts_ch = Channel.fromPath(params.countsTranscripts)

ch_leafcutter_cluster = file(params.leafcutter_cluster)
ch_phenotype_table = file(params.phenotype_table)

/* parameters declared in json file */

phenotype_PCs_sQTL_ch = params.phenotype_PCs_sQTL

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum

/* Genotype common sample */
filter_lst_gene = file( params.sampleIND_gene ).readLines()

filter_lst_transcript = file( params.sampleIND_transcript ).readLines()

filter_lst_splice = file( params.sampleIND_splice ).readLines()

/* calling different process modules */

include { geneTSVfiles } from './modules_dsl2/GeneCountMatrices_TenPercent'

include { mergeNormalizedGeneCountMatrices } from './modules_dsl2/GeneCountMatrices_TenPercent'

include { transcriptGTFfiles } from './modules_dsl2/TranscriptCountMatrices_TenPercent'

include { mergeNormalizedTranscriptCountMatrices } from './modules_dsl2/TranscriptCountMatrices_TenPercent'

include { rnaspliceFilterJunc } from './modules_dsl2/Quantification_SpliceVariants'

/* NOTE: The default PCs for sQTL were defined as 10 but user can modify it in this script */
include { leafcutter_cluster_junctions } from './modules_dsl2/Quantification_SpliceVariants'


workflow {

  /* create a channel to filter RNAseq gene count sample having corresponding genotypes */ 
  gene_St_counts_ch.filter{ it.getName() in filter_lst_gene }
                .set{genefiles_ch}
   
  geneTSVfiles(genefiles_ch)

  mergeNormalizedGeneCountMatrices(geneTSVfiles.out.toList())



  /* create a channel to filter RNAseq transcript count sample having corresponding genotypes */ 
   transcript_counts_ch.filter{ it.getName() in filter_lst_transcript }
                .set{transcriptfiles_ch}

   transcriptGTFfiles(transcriptfiles_ch)

   mergeNormalizedTranscriptCountMatrices(transcriptGTFfiles.out.toList())

  

   /* create a channel to filter junc files if there corresponding .junc sample in found in genotype data "filter_lst" */
   junc_liver_ch.filter{ it.getName() in filter_lst_splice }
                .set{juncfiles_ch}

   rnaspliceFilterJunc(juncfiles_ch)

   leafcutter_cluster_junctions(rnaspliceFilterJunc.out.map{ it.toString() }.collectFile(name: "junction_files.txt", newLine: true),ch_leafcutter_cluster,ch_phenotype_table,intron_length_minimum_ch,intron_length_maximum_ch, phenotype_PCs_sQTL_ch)

  }
