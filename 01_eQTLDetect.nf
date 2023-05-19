nextflow.enable.dsl=2

/*   ------------------- eQTL Nextflow pipeline script01: RNAseq read trimming, alignment and quantification -------------------------- */

/* $projectDir is the current working directory where the script is located */

/* NOTE : Take care of input and output paths before running the script  */


                                    /* --Input: Input RNAseq reads-- */

                                 /* Declare paired-end RNAseq fastq files*/

/*NOTE: Care should be taken to have unique sample ID for each sample before first "_" to extract each sample pair properly  */

/*  1. RNAseq: paired-end library first stranded */
params.reads="$projectDir/Demo_RNAseqData_BovReg/*_R{1,2}_subset.fastq.gz"

/*  2. RNAseq: paired-end library unstranded */
//empty

/*  3. RNAseq: paired-end library second stranded */

// empty

                       /* Declare singel-end RNAseq fastq files*/

/* NOTE: Care should be taken to have unique sample ID for each sample before first "_" to extract each sample properly  */

/*  4. RNAseq: single-end library first stranded */
//empty

/*  5. RNAseq: single-end library unstranded */
//empty

/*  6. RNAseq: single-end library secong stranded */

// empty

                                  /*  --Output: Path to output files-- */

params.outdir = "$projectDir/Output"

                             /*  --Path to supporting files:  gtf, starindex, fasta and software--  */

/* Annotaion file in gtf format*/
params.gtf = "$projectDir/Reference_Genome/Bos_taurus.ARS-UCD1.2.109.gtf"

/* path to STAR reference index created in script 00 based on referece genome */
params.starindex = "$projectDir/GenomeIndex_Ensembl109/star/"

/* trimmomatic file to perform trimming */
params.trimmomaticjar = "$projectDir/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"

/* trimming adapter for pair-end reads */
params.adapter_PE = "$projectDir/Softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

/* trimming adapter for single-end reads */
params.adapter_SE = "$projectDir/Softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"


                                        def inputfiles() {


                                         log.info """\
                                                 R N A S E Q - N F   P I P E L I N E  F O R eQTL DSL2
                                                 ======================================================
                                                 Input        : ${params.reads}
                                                 output       : ${params.outdir}
                                                 gtf          : ${params.gtf}

                                                 """
                                                 .stripIndent()
                                          }

                                        inputfiles()


                                        /** --Channel objects-- **/

/* singel channel objects*/

Channel
       .fromPath(params.starindex)
       .ifEmpty{error "Connot find reference genome: ${params.starindex}"}
       .set {star_index}

trimmomaticjar_ch = file(params.trimmomaticjar)
adapter_PE_ch=file(params.adapter_PE)
adapter_SE_ch=file(params.adapter_SE)


Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set {gtf_ARS}
         

                             /* --Channel for input RNAseq reads-- */

     /* Reads with first stranded library */

read_pairs_ch = Channel.fromFilePairs(params.reads)     



   /* Reads with unstranded library */

//Empty

   /* Reads with second stranded library */

//Empty


/* parameters declared in json file */

anchor_length_ch = params.anchor_length

intron_length_minimum_ch = params.intron_length_minimum

intron_length_maximum_ch = params.intron_length_maximum


                             /*  --Calling modules-- */

include { fastqc_PE } from './modules_dsl2/QualityCheckandQualityControl'

include { fastqc_SE } from './modules_dsl2/QualityCheckandQualityControl'

include { trimmomatic } from './modules_dsl2/QualityCheckandQualityControl'

include { trimmomatic_singleEnd} from './modules_dsl2/QualityCheckandQualityControl'

include { starAlign_GeneTranscript } from './modules_dsl2/Alignment'

include { starAlign_singleEnd_GeneTranscript } from './modules_dsl2/Alignment'

include { starAlign_Splicing } from './modules_dsl2/Alignment'

include { starAlign_singleEnd_Splicing } from './modules_dsl2/Alignment'

include { starAlign_unpaired } from './modules_dsl2/Alignment'

include {samtools_merge_forStringtie} from './modules_dsl2/Indexing'

include {samtools_merge_forleafcutter } from './modules_dsl2/Indexing'

include {samtools_singleEnd_forleafcutter } from './modules_dsl2/Indexing'
                              
include { stringtieQuant_library_firststranded } from './modules_dsl2/Quantification_GeneTranscript'

include { stringtieQuant_library_unstranded } from './modules_dsl2/Quantification_GeneTranscript'

include { stringtieQuant_library_secondstranded } from './modules_dsl2/Quantification_GeneTranscript'

include { leafcutter_bamtojunc_library_firststranded } from './modules_dsl2/Quantification_SpliceVariants'

include { leafcutter_bamtojunc_library_unstranded } from './modules_dsl2/Quantification_SpliceVariants'

include { leafcutter_bamtojunc_library_secondstranded } from './modules_dsl2/Quantification_SpliceVariants'


/* Calling different sub-workflows in the main workflow 
   
    NOTE: Users can run only one or any sub-workflow combinations 
    by muting unwanted sub-worflows with comments // at the beginnning of workflow declaration 
    for example: //SINGELEND_READS_FIRST_STRANDED_LIBRARY(....) */



workflow {

    main:

     PAIREDEND_READS_FIRST_STRANDED_LIBRARY(read_pairs_ch, trimmomaticjar_ch,adapter_PE_ch,gtf_ARS,star_index,anchor_length_ch,intron_length_minimum_ch,intron_length_maximum_ch)

     //SINGELEND_READS_FIRST_STRANDED_LIBRARY(sample_ch,trimmomaticjar_ch,adapter_SE_ch,gtf_ARS,star_index,anchor_length_ch,intron_length_minimum_ch,intron_length_maximum_ch)

     //PAIREDEND_READS_UNSTRANDED_LIBRARY(sample_ch,trimmomaticjar_ch,adapter_PE_ch,gtf_ARS,star_index,anchor_length_ch,intron_length_minimum_ch,intron_length_maximum_ch )

     //SINGELEND_READS_UNSTRANDED_LIBRARY(sample_ch,trimmomaticjar_ch,adapter_SE_ch,gtf_ARS,star_index,anchor_length_ch,intron_length_minimum_ch,intron_length_maximum_ch )

     //PAIREDEND_READS_SECOND_STRANDED_LIBRARY(sample_ch, trimmomaticjar_ch,adapter_PE_ch,gtf_ARS,star_index,anchor_length_ch,intron_length_minimum_ch,intron_length_maximum_ch)

     //SINGELEND_READS_SECOND_STRANDED_LIBRARY(sample_ch,trimmomaticjar_ch,adapter_SE_ch,gtf_ARS,star_index,anchor_length_ch,intron_length_minimum_ch,intron_length_maximum_ch)   
    
}



/* This workflow takes the input RNAseq reads which are paired-end and having stranded library fr-firststrand  */

workflow PAIREDEND_READS_FIRST_STRANDED_LIBRARY {

        take:

          paired_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

        main:

         fastqc_PE(paired_reads)

         trimmomatic(paired_reads,trimmjar_file,adapter_file)

         starAlign_GeneTranscript(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

         starAlign_Splicing(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

         starAlign_unpaired(trimmomatic.out.forward_unpaired_trim.join(trimmomatic.out.revese_unpaired_trim),gtf_file.collect(),star_index.collect())
   
         samtools_merge_forStringtie(starAlign_GeneTranscript.out.tomerge_paired_ch.join(starAlign_unpaired.out.tomerge_unpairedR1_ch).join(starAlign_unpaired.out.tomerge_unpairedR2_ch))

         samtools_merge_forleafcutter(starAlign_Splicing.out.tomerge_paired_Leaf_ch.join(starAlign_unpaired.out.tomerge_unpairedLeafR1_ch).join(starAlign_unpaired.out.tomerge_unpairedLeafR2_ch))
        
         stringtieQuant_library_firststranded(samtools_merge_forStringtie.out.tostringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_firststranded(samtools_merge_forleafcutter.out.toleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum)

}


/* This workflow takes the input RNAseq reads which are single-end and having stranded library fr-firststrand  */

workflow SINGELEND_READS_FIRST_STRANDED_LIBRARY {

        take:

          singelend_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

        main:

         fastqc_SE(singelend_reads)

         trimmomatic_singleEnd(singelend_reads,trimmjar_file,adapter_file)

         starAlign_singleEnd_GeneTranscript(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         starAlign_singleEnd_Splicing(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         samtools_singleEnd_forleafcutter(starAlign_singleEnd_Splicing.out.samindex_ch)
    
         stringtieQuant_library_firststranded(starAlign_singleEnd_GeneTranscript.out.stringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_firststranded(samtools_singleEnd_forleafcutter.out.tosingelEndleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum) 

}

/* This workflow takes the input RNAseq reads which are paired-end and having unstranded library */

workflow PAIREDEND_READS_UNSTRANDED_LIBRARY {

        take:

          paired_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

        main:

         fastqc_PE(paired_reads)

         trimmomatic(paired_reads,trimmjar_file,adapter_file)

         starAlign_GeneTranscript(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

         starAlign_Splicing(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

         starAlign_unpaired(trimmomatic.out.forward_unpaired_trim.join(trimmomatic.out.revese_unpaired_trim),gtf_file.collect(),star_index.collect())
   
         samtools_merge_forStringtie(starAlign_GeneTranscript.out.tomerge_paired_ch.join(starAlign_unpaired.out.tomerge_unpairedR1_ch).join(starAlign_unpaired.out.tomerge_unpairedR2_ch))

         samtools_merge_forleafcutter(starAlign_Splicing.out.tomerge_paired_Leaf_ch.join(starAlign_unpaired.out.tomerge_unpairedLeafR1_ch).join(starAlign_unpaired.out.tomerge_unpairedLeafR2_ch))
        
         stringtieQuant_library_unstranded(samtools_merge_forStringtie.out.tostringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_unstranded(samtools_merge_forleafcutter.out.toleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum)  
}


/* This workflow takes the input RNAseq reads which are single-end and having unstranded library */

workflow SINGELEND_READS_UNSTRANDED_LIBRARY {

        take:

          singelend_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

        main:

         fastqc_SE(singelend_reads)

         trimmomatic_singleEnd(singelend_reads,trimmjar_file,adapter_file)

         starAlign_singleEnd_GeneTranscript(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         starAlign_singleEnd_Splicing(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         samtools_singleEnd_forleafcutter(starAlign_singleEnd_Splicing.out.samindex_ch)
    
         stringtieQuant_library_unstranded(starAlign_singleEnd_GeneTranscript.out.stringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_unstranded(samtools_singleEnd_forleafcutter.out.tosingelEndleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum) 
}


/* This workflow takes the input RNAseq reads which are paired-end and having stranded library fr-secondstrand */

workflow PAIREDEND_READS_SECOND_STRANDED_LIBRARY {

        take:

          paired_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

        main:

         fastqc_PE(paired_reads)

         trimmomatic(paired_reads,trimmjar_file,adapter_file)

         starAlign_GeneTranscript(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

         starAlign_Splicing(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

         starAlign_unpaired(trimmomatic.out.forward_unpaired_trim.join(trimmomatic.out.revese_unpaired_trim),gtf_file.collect(),star_index.collect())
   
         samtools_merge_forStringtie(starAlign_GeneTranscript.out.tomerge_paired_ch.join(starAlign_unpaired.out.tomerge_unpairedR1_ch).join(starAlign_unpaired.out.tomerge_unpairedR2_ch))

         samtools_merge_forleafcutter(starAlign_Splicing.out.tomerge_paired_Leaf_ch.join(starAlign_unpaired.out.tomerge_unpairedLeafR1_ch).join(starAlign_unpaired.out.tomerge_unpairedLeafR2_ch))
        
         stringtieQuant_library_secondstranded(samtools_merge_forStringtie.out.tostringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_secondstranded(samtools_merge_forleafcutter.out.toleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum)  
}


/* This workflow takes the input RNAseq reads which are single-end and having stranded library fr-secondstrand */

workflow SINGELEND_READS_SECOND_STRANDED_LIBRARY {

        take:

          singelend_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

        main:

         fastqc_SE(singelend_reads)

         trimmomatic_singleEnd(singelend_reads,trimmjar_file,adapter_file)

         starAlign_singleEnd_GeneTranscript(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         starAlign_singleEnd_Splicing(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         samtools_singleEnd_forleafcutter(starAlign_singleEnd_Splicing.out.samindex_ch)
    
         stringtieQuant_library_secondstranded(starAlign_singleEnd_GeneTranscript.out.stringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_secondstranded(samtools_singleEnd_forleafcutter.out.tosingelEndleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum) 

}










                                              /* ****END****** */




