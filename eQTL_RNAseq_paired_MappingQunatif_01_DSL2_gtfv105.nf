nextflow.enable.dsl=2

/* $projectDir is the current working directory /disk2/chitneedi/NextFlow */

/**           ****************** NOTE ******************                           */
/* This script is designed to run the paired-end input RNA-Seq reads              */
/* Take care of input and output paths before running the script                  */
/* This script performs Quality check, trimming, Alignment and quantification at gene level, transcript level and splicing regions */

/*  Script parameters */

                                       /*  FBN  */  
                        /* NOTE: Take care of StringTie strandedness */
/* Liver  */                                       
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Liver/*R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Liver"

/* Rumen  */  
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Rumen/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Rumen"

/* Muscle  */ 
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Muscle/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Muscle"

/* Jejunum  */ 
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Jejunum/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Jejunum"

/* MammaryGland  */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/MammaryGland/HL/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/MammaryGland/HL"

//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/MammaryGland/VL/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/MammaryGland/VL"

/* Blood  */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Blood/00h/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Blood/00h"

//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Blood/24h/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Blood/24h"

//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN/Blood/96h/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Blood/96h"

                                     /*  INRAE  */ 
                        /* NOTE: Take care of StringTie strandedness */
/* Muscle */ 
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/INRAE/Merged_INRAE_RNAseq/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/INRAE_eQTL/INRAE_RNASeq"

/* Heart */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/INRAE/Merged_INRAE_RNASeq_Heart/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/INRAE_eQTL/INRAE_RNASeq_Heart"

/* Liver */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/INRAE/Merged_INRAE_RNASeq_Liver/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/INRAE_eQTL/INRAE_RNASeq_Liver"

/* Cheek-Muscle" */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/INRAE/Merged_INRAE_RNASeq_HOL_Muscle-cheek/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/INRAE_eQTL/INRAE_RNAseq_HOL_Muscle-cheek"


                                    /*  LUKE  */ 
                        /* NOTE: Take care of StringTie strandedness */
   /* Adipose */                     
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/LUKE/bovreg/Bovine_adipose/*_{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_adipose"
  /* Adipose */ 
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/LUKE/bovreg/Bovine_adipose2/*_{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_adipose2"

  /* Liver */ 
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/LUKE/bovreg/Bovine_liver/*_{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_liver"
 /* Liver */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/LUKE/bovreg/Bovine_liver2/*_{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_liver2"
 /* Mammary gland */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/LUKE/bovreg/Bovine_mammaryGland/*_{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_mammaryGland"
/* Bovine papilli */
//params.reads="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/LUKE/bovreg/Bovine_papilli/*_{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_papilli"


                                          /*  UAL  */  
                        /* NOTE: Take care of StringTie strandedness */

                             /* Graham_Plastow_data */
                       /* NOTE: Take care of StringTie strandedness */
//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/UAL/Graham_Plastow_BRD_RNA-Seq_data/MCG_4362_20170922/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/UAL_eQTLGraham_Plastow_BRD_RNA-Seq_data/MCG_4362_20170922"

//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/UAL/Graham_Plastow_BRD_RNA-Seq_data/raw_fastq1/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/UAL_eQTL/Graham_Plastow_BRD_RNA-Seq_data/raw_fastq1"

                                      
                                       /*  FBN_add_data_sets   */

//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN_add_data_sets/Clust_Jejunum_Merged/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL_add_data_sets/Clust_Jejunum"

//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN_add_data_sets/Clust_Rumen_Merged/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL_add_data_sets/Clust_Rumen"

//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/FBN_add_data_sets/Koch_Blood_Merged/*_R{1,2}.fastq.gz"
//params.outdir ="$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL_add_data_sets/Koch_Blood"
                   
                                   /*Ag Vic data sets */

//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/AgVic_Preprocessing/*_blood_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/AgVic_eQTL/Blood" 

//params.reads = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/Data/AgVic_Preprocessing/2*_milk_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/AgVic_eQTL/Milk"




/* - The Ensembl gtf version 105 was used to check the differences in gene and transcript counts with old Alignment
   - The v105 will only used from stringtie quantification and subsequent steps.

  */

params.gtf = "$projectDir/BovReg_Fasta/Bos_taurus.ARS-UCD1.2.105.gtf"


/* Liver Aligned reads as input for quantifiaction */

params.FBN_paired = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Liver/star_aligned/*_paired_Aligned.sortedByCoord.out.bam"
params.FBN_unpairedR1 = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Liver/star_aligned/*_unpaired_R1_Aligned.sortedByCoord.out.bam"
params.FBN_unpairedR2 = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/FBN_eQTL/Liver/star_aligned/*_unpaired_R2_Aligned.sortedByCoord.out.bam"


params.GIGA_singelEnd = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/GIGA_eQTL/Liver/star_aligned/*_Aligned.sortedByCoord.out.bam"

params.LUKE_paired = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_liver*/star_aligned/*_paired_Aligned.sortedByCoord.out.bam"
params.LUKE_unpairedR1 = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_liver*/star_aligned/*_unpaired_R1_Aligned.sortedByCoord.out.bam"
params.LUKE_unpairedR2 = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/LUKE/Bovine_liver*/star_aligned/*_unpaired_R2_Aligned.sortedByCoord.out.bam"

params.UAL_singelEnd = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/UAL_eQTL/Changxi_Li_RNA/LIVER_Tissue/star_aligned/*_Aligned.sortedByCoord.out.bam"


params.outdir = "$projectDir/../winpro/I3-BovReg_Genombiologie/WP4_eQTL/eQTL_analyses/Liver_eQTL/GTF_v105"

//params.gtf = "$projectDir/BovReg_Fasta/BovReg_mRNA_totalRNA_CAGE.gff"
params.starindex = "$projectDir/GenomeIndex_BovReg/star/"
params.fasta = "$projectDir/BovReg_Fasta/ARS-UCD1.2_Btau5.0.1Y.fa"
params.outgff = "$projectDir/ARS-UCD_Fasta"
params.phenotype_table = "$projectDir/leafcutter/scripts/prepare_phenotype_table.py"
params.leafcutter_cluster = "$projectDir/leafcutter/scripts/leafcutter_cluster_Bovine_regtools.py"
//params.dexseq = "$projectDir/dexseq_Dockerfile/dexseq_prepare_annotation2.py"
params.trimmomaticjar = "$projectDir/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"
params.adapter = "$projectDir/Softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"



def inputfiles() {


 log.info """\
         R N A S E Q - N F   P I P E L I N E  F O R eQ T L DSL2
         ======================================================
         test          : ${params.FBN_paired}
         output       : ${params.outdir}
         fasta_used   : ${params.fasta}
         gtf          : ${params.gtf}

         """
         .stripIndent()
  }

inputfiles()


/* Channel objects */

/* single channel objects*/
/*
Channel
       .fromPath(params.starindex)
       .ifEmpty{error "Connot find reference genome: ${params.starindex}"}
        .set {star_index}
*/
    

trimmomaticjar_ch = file(params.trimmomaticjar)
adapter_ch=file(params.adapter)


ch_leafcutter_cluster = file(params.leafcutter_cluster)
ch_phenotype_table = file(params.phenotype_table)

//read_pairs_ch = Channel.fromFilePairs(params.reads)


Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set {gtf_ARS}
         

/* Channel for aligned reads */
     /* Qunatification with --rf */

Channel
        .fromPath(params.FBN_paired)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_paired_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {paired_FBN}
       

Channel
        .fromPath(params.FBN_unpairedR1)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_unpaired_R1_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {unpairedR1_FBN}

Channel
        .fromPath(params.FBN_unpairedR2)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_unpaired_R2_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {unpairedR2_FBN}

Channel
        .fromPath(params.GIGA_singelEnd)
         .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {singelEnd_GIGA}


       /* Qunatification without  --rf */

Channel
        .fromPath(params.LUKE_paired)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_paired_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {paired_LUKE}
       


Channel
        .fromPath(params.LUKE_unpairedR1)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_unpaired_R1_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {unpairedR1_LUKE}
       
 
Channel
        .fromPath(params.LUKE_unpairedR2)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_unpaired_R2_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {unpairedR2_LUKE}
          




Channel
        .fromPath(params.UAL_singelEnd)
        .map { path ->
             def content = path 
             def name = path.baseName.replaceAll('_Aligned.sortedByCoord.out','')
             [ name, content]  }
        .set {singelEnd_UAL}



// Temp channel objects

//tomerge_ch = Channel.fromPath(params.samtools_input)
       
//stringtie_ch = Channel.fromPath(params.stringtie_input)
       

/* calling different process modules */

include { fastqc } from './modules_dsl2/QualityCheckandQualityControl'

include { trimmomatic } from './modules_dsl2/QualityCheckandQualityControl'

include { trimmomatic_singleEnd} from './modules_dsl2/QualityCheckandQualityControl'

include { starAlign_GeneTranscript } from './modules_dsl2/Alignment'

include { starAlign_Splicing } from './modules_dsl2/Alignment'

include { starAlign_unpaired } from './modules_dsl2/Alignment'

include {samtools_merge_forStringtie_with_rf} from './modules_dsl2/Indexing'

include {samtools_merge_forStringtie_without_rf} from './modules_dsl2/Indexing'

include {samtools_merge_forleafcutter } from './modules_dsl2/Indexing'
                              
include { stringtieQuant_with_rf } from './modules_dsl2/Quantification_GeneTranscript'

include { stringtieQuant_without_rf } from './modules_dsl2/Quantification_GeneTranscript'

include { leafcutter_junctions } from './modules_dsl2/Quantification_SpliceVariants'





workflow {

    main:
     //fastqc(read_pairs_ch)

    //trimmomatic(read_pairs_ch,trimmomaticjar_ch,adapter_ch)

    //starAlign_GeneTranscript(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_ARS.collect(),star_index.collect())

    //starAlign_Splicing(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_ARS.collect(),star_index.collect())

    //starAlign_unpaired(trimmomatic.out.forward_unpaired_trim.join(trimmomatic.out.revese_unpaired_trim),gtf_ARS.collect(),star_index.collect())
  

    fbn_ch = paired_FBN.join(unpairedR1_FBN).join(unpairedR2_FBN)
    samtools_merge_forStringtie_with_rf(fbn_ch)

    luke_ch = paired_LUKE.join(unpairedR1_LUKE).join(unpairedR2_LUKE) 
    samtools_merge_forStringtie_without_rf(luke_ch)


    //samtools_merge_forleafcutter(starAlign_Splicing.out.tomerge_paired_Leaf_ch.join(starAlign_unpaired.out.tomerge_unpairedLeafR1_ch).join(starAlign_unpaired.out.tomerge_unpairedLeafR2_ch))
   

   
    liver_with_rf = samtools_merge_forStringtie_with_rf.out.tostringtie_with_rf.concat(singelEnd_GIGA)
   
    
    stringtieQuant_with_rf(liver_with_rf,gtf_ARS.collect())


    liver_without_rf = samtools_merge_forStringtie_without_rf.out.tostringtie_without_rf.concat(singelEnd_UAL)
    

    stringtieQuant_without_rf(liver_without_rf,gtf_ARS.collect())





    //leafcutter_junctions(samtools_merge_forleafcutter.out)

    //leafcutter_cluster_junctions(leafcutter_junctions.out.map { it.toString() }.collectFile(name: "junction_files.txt", newLine: true),ch_leafcutter_cluster,ch_phenotype_table)
    
}














                                       /*--------unused processes----------*/


                    /*
                     * STEP 7 - Quantify exon expression - featureCounts (read count and exon counts )
                     */

/*
Kerimov et al.  "we first created exon annotation file (GFF) using GENCODE V30 reference
transcriptome annotations and dexseq_prepare_annotation.py script from the DEXSeq [19]
package. We then used the aligned RNA-seq BAM files from the gene expression quantification
and featureCounts with flags ‘-p -t exonic_part -s ${direction} -f -O’ to count the
number of reads overlapping each exon".
*/

/*
process makeDexSeqExonGFF {
        publishDir "${params.outgff}", mode: 'copy'
        container 'praveen/pythonhtseq'

        input:
        file gtf from gtf_dexseq
        file dexseq from ch_dexseq

        output:
        file "Bos_taurus.ARS-UCD1.2patchedDEXSeq.gff" into ch_dexseq_gff_count_exons

        file "Bos_taurus.ARS-UCD1.2patchedDEXSeq.gtf" into ch_dexseq_gtf_count_exons

        script:
        //cat $gtf | sed 's/chrM/chrMT/;s/chr//' > ${gtf.baseName}.patched_contigs.gtf
        """
         python $dexseq -f Bos_taurus.ARS-UCD1.2patchedDEXSeq.gtf $gtf Bos_taurus.ARS-UCD1.2patchedDEXSeq.gff
        """
}
*/



/*
process featureCounts_quant {
     container 'genomicpariscentre/featurecounts'
     publishDir "${params.outdir}/featurecountQuantify", mode:'copy'
     input:
      file  bam  from featureCounts_ch
      file gtf from gtf_featureCounts.collect()
      file gff from ch_dexseq_gff_count_exons.collect()
    output:

     file "*" into featurecounts_output_ch

    script:
     samplename = bam[0].toString() - ~/(\.sortedByCoord)?(\.out)?(\.bam)?$/

    //Get read counts
    // the option -f (f specified, read summarization will be performed at featurelevel  (eg.   exon  level).
    // Otherwise,  it  is  performed  at  meta-feature level (eg.  gene level).) has to checked for exon count !!

    """
    featureCounts -T 10 --donotsort -p -s 2 -a $gtf -o ${samplename}_featureCounts_readCount.txt $bam
    featureCounts -T 10 -p -t exonic_part -f -O -s 2 -a $gff -o ${samplename}_featureCounts_exonCount.txt $bam

    """
    //awk '!visitedline[\$0]++ ' ${samplename}_featureCounts_exonCount.txt > ${samplename}_featureCounts_exonCount_Mod.txt
    //awk -F'\t' -v OFS='\t' '{ for (i = 1 ; i <= NF ; i++) gsub(/;.* /,"",\$i)}1' ${samplename}_featureCounts_exonCount.txt > ${samplename}_featureCounts_exonCount_Mod.txt

}
*/