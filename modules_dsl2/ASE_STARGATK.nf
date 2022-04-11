nextflow.enable.dsl=2
/*   ------------------- eQTL Nextflow pipeline script01: RNA mapping and read quantification -------------------------- */


/*  Script parameters */

//params.reads ="$baseDir/mntDirPrav/disk2mnt/liver_FBN/B01*_liver_R{1,2}.fastq.gz"
params.fasta = "$baseDir/BovReg_Fasta/ARS-UCD1.2_Btau5.0.1Y.fa"
params.fai = "$baseDir/BovReg_Fasta/ARS-UCD1.2_Btau5.0.1Y.fa.fai"
params.dict = "$baseDir/BovReg_Fasta/ARS-UCD1.2_Btau5.0.1Y.dict"
params.gtf = "$baseDir/BovReg_Fasta/BovReg_mRNA_totalRNA_CAGE.gff"
params.starindex = "$baseDir/GenomeIndex_BovReg/star/"
params.outdir = "$baseDir/../../../disk2/ASE_Data/Results/"
//params.outdirtmp = "$baseDir/mntDirPrav/disk2mnt/Inter_results_eQTL_MappingQunatif"
//----The Vcf file downloaded from https://ftp.ensembl.org/pub/release-102/variation/vcf/bos_taurus/
// --- The "indels" were removed from the vcf file was then recoded and indexed
// Command: ../../softwareTools/vcftools_0.1.13/bin/vcftools --gzvcf bos_taurus_102.vcf.gz  --remove-indels -c --recode | bgzip -c > bos_taurus_102.filteredIndels.recode.vcf.gz
//params.inputvcf = "$baseDir/eQTL_Genotype_vcf/liver_Test_Genotypes/ImpuHDktoWGS_RNA_liver_Chr_25_UpIds.vcf"
params.inputvcf = "$baseDir/../../../disk2/ASE_Data/Genotype_FBN/ImpuHDktoWGS_segfam_ChromMas_NUDA_Ars1.2_NdDidierfile_1043IND_Chr11.vcf.gz"
params.outgff = "$baseDir/ARS-UCD_Fasta"
params.phenotype_table = "$baseDir/leafcutter/scripts/prepare_phenotype_table.py"
params.leafcutter_cluster = "$baseDir/leafcutter/scripts/leafcutter_cluster_Bovine_regtools.py"
params.dexseq = "$baseDir/dexseq_Dockerfile/dexseq_prepare_annotation2.py"
params.trimmomaticjar = "$baseDir/Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"
params.adapter = "$baseDir/Softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
params.inputreads = "$baseDir/../../../disk2/ASE_Data/RNAseqData/H11_HL_R{1,2}.clean.fastq.gz"
params.picardjar = "$baseDir/Softwares/Picard-2.25.5/picard.jar"
//params.trimmomaticRead = "$baseDir/results_eQTL_MappingQunatif/Trimmoatic_output/*_liver_R{1,2}.clean.fastq.gz"

log.info """\
         R N A S E Q - N F   P I P E L I N E  F O R eQ T L
         =================================================
         gtf          : ${params.gtf}
         reads        : ${params.inputreads}
         fasta        : ${params.fasta}

         """
         .stripIndent()


/* Channel objects */

/* single channel objects*/

Channel
       .fromPath(params.starindex)
       .ifEmpty{error "Connot find reference genome: ${params.starindex}"}
        .set {star_index}
trimmomaticjar_ch = file(params.trimmomaticjar)
adapter_ch=file(params.adapter)




Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .set {gtf_ARS}

ch_leafcutter_cluster = file(params.leafcutter_cluster)
ch_phenotype_table = file(params.phenotype_table)
picardjar_ch = file(params.picardjar)

//ch_dexseq = file(params.dexseq)

/* multi channel objects*/
/*
Channel
      .fromFilePairs(params.trimmomaticRead)
      .ifEmpty { error "Cannot find any reads matching: ${params.trimmomaticRead}" }
      .set {read_splicevar_ch}
      */

Channel
    .fromFilePairs(params.inputreads )
    .ifEmpty { error "Cannot find any reads matching: ${params.inputreads}" }
    .set { inputreads_ch }


Channel
       .fromPath(params.inputvcf)
       .ifEmpty {exit 1, "vcf file not found: ${params.inputvcf}"}
       .set {inputvcf_ch}

Channel
      .fromPath(params.fasta)
      .ifEmpty { error "Cannot find any reads matching: ${params.fasta}" }
      .set {fasta_ch}

Channel
      .fromPath(params.fai)
      .ifEmpty { error "Cannot find any reads matching: ${params.fai}" }
      .set {fai_ch}

Channel
      .fromPath(params.dict)
      .ifEmpty { error "Cannot find any reads matching: ${params.dict}" }
      .set {dict_ch}





/* Join the cleaned forward and reverse strand read pairs for STAR alignment*/
 //strands_ch = forward_strand_ch.merge(reverse_strand_ch) 


/*
// Trimmomatic report extracted using the following command
//grep -e "Input Read Pairs:" $baseDir/results_eQTL_MappingQunatif/Trimmoatic_output/*_liver_trimomatic.log | awk '{print $1,$4,$7$8,$12$13,$17$18,$20$21}' | sed 's/_trimomatic.log:Input//g' | sed '1iSample Input_Read_Pairs Both_Surviving Forward_Only_Surviving Reverse_Only_Surviving Dropped' | tr ' ' '\t'  >> Trimmomatic_output_report.txt
*/
                        /*
                         * STEP 2 - Read alignment with STAR and the output BAM files are filtered for wasp filter 
                                    */
/*
 process starAlign {
    container 'mgibio/star'
    publishDir "${params.outdir}/star_aligned_ASE_WASP", mode:'copy'

    input:
    //set sample_id, file(reads) from read_pairs2_ch
    set file(read1), file(read2)from forward_strand_ch.merge(reverse_strand_ch)
    file gtf from gtf_star.collect()
    file index from star_index.collect()
    file vcf from inputvcf.collect()
    //set sample_id, file (read_splicevar) from read_splicevar_ch
    output:
    file("*Log.final.out") into star_aligned1
    file "${prefix}_NoWASP_Aligned.sortedByCoord.out.bam" into samwasp_ch
    //file "${sample_id}_Aligned.toTranscriptome.out.bam" into star_Transc_bam
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    //file "${sample_id}_leafcutter_Aligned.out.bam" into leafcutter_bam_ch

    script:
    prefix = read1[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/



    """
      STAR --genomeDir $index \
             --readFilesIn ${read1} ${read2} \
             --runThreadN 10 \
             --outSAMtype BAM SortedByCoordinate \
             --readFilesCommand zcat \
             --outFileNamePrefix ${prefix}_NoWASP_
    """
     
  WASP filter
     STAR --genomeDir $index \
             --readFilesIn ${read1} ${read2} \
             --runThreadN 10 \
             --outSAMtype BAM SortedByCoordinate \
             --waspOutputMode SAMtag \
             --varVCFfile ${vcf}\
             --readFilesCommand zcat \
             --outFileNamePrefix ${prefix}_

     STAR --genomeDir $index \
             --sjdbGTFfile $gtf \
             --readFilesIn  $read_splicevar \
             --twopassMode Basic \
             --outSAMstrandField intronMotif \
             --readFilesCommand zcat \
             --outSAMtype BAM Unsorted \
             --outFileNamePrefix ${sample_id}_leafcutter_ \
             --runThreadN 10
        """
        
}
*/
                       /*
                         * STEP 3 - Filter WASp filter tags  
                         */
/*
process samtoolsFilterWASPtags {
    container 'biocontainers/samtools:v1.9-4-deb_cv1'
    publishDir "${params.outdir}/leafcutter_bamtojunc", mode: 'copy'

    input:
    file bamwasp from samwasp_ch
    output:

    file "${samplename}_WASPFiltered.bam" into bam_rg_readCount_ch
   
    script:
 

  samplename = bamwasp[0].toString() - ~/(_Aligned)?(\.sortedByCoord)?(\.out)?(\.bam)?$/

  """
   samtools view -H ${bamwasp} >  ${samplename}_header.sam

   samtools view ${bamwasp} | grep  'vW:i:1' | grep -w 'NH:i:1' > ${samplename}_noheader.sam

   cat  ${samplename}_header.sam ${samplename}_noheader.sam | samtools view -Sb  > ${samplename}_WASPFiltered.bam

   rm ${samplename}_header.sam ${samplename}_noheader.sam

  """
}


process Addorreplacereadgroup{

  input:
        file picardjar from picardjar_RG_ch

        file input_bam from inputreads_ch

 output:
        file "${prefix}_RG.bam" into  bam_rg_readCount_ch

 script:
      
      prefix = input_bam[0].toString() - ~/(_WASPFiltered)?(\.bam)?$/ 

      """
      java -jar $picardjar AddOrReplaceReadGroups \
       I=${input_bam} \
       O=${prefix}_RG.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20
      """
}
*/


/* calling different process modules */

include { starAlign1 } from './modules_dsl2/Alignment_ASE'

include { starAlign2 } from './modules_dsl2/Alignment_ASE'

include { addorreplacereadgroup } from './modules_dsl2/ASE_RNAFilters'

include { markDuplicates } from './modules_dsl2/ASE_RNAFilters'

include { splitNCigarReads } from './modules_dsl2/ASE_RNAFilters'

include { indexvcf_ASE } from './modules_dsl2/IndexVcf'

include { readCounter } from './modules_dsl2/ASEReadCounter'

workflow {

    main:
    
    starAlign1(inputreads_ch,star_index.collect(),gtf_ARS.collect())

    starAlign2(inputreads_ch,star_index.collect(),gtf_ARS.collect())

    addorreplacereadgroup(picardjar_ch,starAlign1.out.ase_bam_ch1, starAlign2.out.ase_bam_ch2 )

    markDuplicates(picardjar_ch, addorreplacereadgroup.out.bam_rg_ch1, addorreplacereadgroup.out.bam_rg_ch2)
    
    // Adding spliNCigarReads dosen't make much difference and this step is alos time consuming 
    //splitNCigarReads(fasta_ch.collect(),fai_ch.collect(),dict_ch.collect(), markDuplicates.out.markDup_ch1,markDuplicates.out.markDup_ch2)

    indexvcf_ASE(inputvcf_ch)

    //readCounter(fasta_ch.collect(),fai_ch.collect(),dict_ch.collect(), indexvcf_ASE.out.vcf.join(indexvcf_ASE.out.vcf_index), splitNCigarReads.out.cigarReads_ch1, splitNCigarReads.out.cigarReads_ch2)
    
   readCounter(fasta_ch.collect(),fai_ch.collect(),dict_ch.collect(), indexvcf_ASE.out.vcf.collect().join(indexvcf_ASE.out.vcf_index.collect()),markDuplicates.out.markDup_ch1, markDuplicates.out.markDup_ch2)

}


                            /* ****END****** */
