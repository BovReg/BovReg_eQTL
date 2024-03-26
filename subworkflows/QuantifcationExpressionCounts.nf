/*
========================================================================================
    Include Modules
========================================================================================
*/

include { stringtieQuant_library_firststranded } from '../modules_dsl2/Quantification_GeneTranscript'

include { stringtieQuant_library_secondstranded} from '../modules_dsl2/Quantification_GeneTranscript'

include { stringtieQuant_library_unstranded} from '../modules_dsl2/Quantification_GeneTranscript'

include { leafcutter_bamtojunc_library_firststranded } from '../modules_dsl2/Quantification_SpliceVariants'

include { leafcutter_bamtojunc_library_secondstranded } from '../modules_dsl2/Quantification_SpliceVariants'

include { leafcutter_bamtojunc_library_unstranded } from '../modules_dsl2/Quantification_SpliceVariants'

include {tissuewise_extractGenotype} from '../modules_dsl2/Tissuewise_Extract_Merge_Genotypes'

include { geneTSVfiles } from '../modules_dsl2/GeneCountMatrices_TenPercent'

include { mergeNormalizedGeneCountMatrices } from '../modules_dsl2/GeneCountMatrices_TenPercent'

include { transcriptGTFfiles } from '../modules_dsl2/TranscriptCountMatrices_TenPercent'

include { mergeNormalizedTranscriptCountMatrices } from '../modules_dsl2/TranscriptCountMatrices_TenPercent'

include { rnaspliceFilterJunc } from '../modules_dsl2/Quantification_SpliceVariants'

include { filter_Splicecounts } from '../modules_dsl2/Quantification_SpliceVariants'

include { leafcutter_cluster_junctions } from '../modules_dsl2/Quantification_SpliceVariants'



params.firstStranded = false

params.secondStranded = false

params.unStranded = false


workflow QUANTIFICATION_AND_MERGE_EXPRESSION_COUNTS {


      take:
           stringtie_ip

           leafcutter_bamtojunc_ip

           gtf_file

           anchor_length 

           intron_length_minimum

           intron_length_maximum

           genotype_input

           sample_info

          phenotype_PCs_sQTL

          leafcutter_cluster_py

          leafcutter_table_py


      main:


       if ( params.firstStranded ) {

	      stringtieQuant_library_firststranded(stringtie_ip,gtf_file.collect())

        leafcutter_bamtojunc_library_firststranded(leafcutter_bamtojunc_ip,anchor_length,intron_length_minimum,intron_length_maximum)
        
        gene_St_counts           = stringtieQuant_library_firststranded.out.gene_St_counts_ch

        splice_junc_sample_tuple = leafcutter_bamtojunc_library_firststranded.out.junc_ch

        transcript_counts        = stringtieQuant_library_firststranded.out.transcript_counts_ch


        }

       if ( params.secondStranded ) {

        stringtieQuant_library_secondstranded(stringtie_ip,gtf_file.collect())

        leafcutter_bamtojunc_library_secondstranded(leafcutter_bamtojunc_ip,anchor_length,intron_length_minimum,intron_length_maximum)
        
        gene_St_counts           = stringtieQuant_library_secondstranded.out.gene_St_counts_ch

        splice_junc_sample_tuple = leafcutter_bamtojunc_library_secondstranded.out.junc_ch

        transcript_counts        = stringtieQuant_library_secondstranded.out.transcript_counts_ch

        

        }

      
       if ( params.unStranded ) {

        stringtieQuant_library_unstranded(stringtie_ip,gtf_file.collect())

        leafcutter_bamtojunc_library_unstranded(leafcutter_bamtojunc_ip,anchor_length,intron_length_minimum,intron_length_maximum)

        gene_St_counts           = stringtieQuant_library_unstranded.out.gene_St_counts_ch

        splice_junc_sample_tuple = leafcutter_bamtojunc_library_unstranded.out.junc_ch  

        transcript_counts        = stringtieQuant_library_unstranded.out.transcript_counts_ch    

        }


         tissuewise_extractGenotype(genotype_input,sample_info.collect())
          
         /* Extract the list of samples belongs to genotype data with RNAseq sample ids*/  
         geno_sample_list=tissuewise_extractGenotype.out.lst_Ind_ch.map{it[1].readLines()}.flatMap()
       
         /* create a channel to filter RNAseq gene count sample having corresponding genotypes */ 
          
         gene_St_counts.join(geno_sample_list).set{genefiles_ch}
          
         geneTSVfiles(genefiles_ch)

         mergeNormalizedGeneCountMatrices(geneTSVfiles.out.toList())

          /* create a channel to filter RNAseq transcript count sample having corresponding genotypes */ 
          
         //transcript_sample_tuple=stringtieQuant_library_firststranded.out.transcript_counts_ch
    
         transcript_counts.join(geno_sample_list).set{transcriptfiles_ch}

         transcriptGTFfiles(transcriptfiles_ch)

         mergeNormalizedTranscriptCountMatrices(transcriptGTFfiles.out.toList())


         /* create a channel to filter junc files if there corresponding .junc sample in found in genotype data */
           
  
         splice_junc_sample_tuple.join(geno_sample_list).set{juncfiles_ch}

         leafcutter_cluster_junctions(juncfiles_ch.map{ it[1].toString() }.collectFile(name: "junction_files.txt", newLine: true),leafcutter_cluster_py,leafcutter_table_py,intron_length_minimum,intron_length_maximum, phenotype_PCs_sQTL)
        
         leafcutter_cluster_junctions.out.spliceCounts
        
         norm_Tr_count =  mergeNormalizedTranscriptCountMatrices.out
        
         norm_Ge_count = mergeNormalizedGeneCountMatrices.out

         splice_Count    = leafcutter_cluster_junctions.out.spliceCounts

         splice_pc    = leafcutter_cluster_junctions.out.splicepcs


         emit:
    
         norm_Tr_count 

         norm_Ge_count 

         splice_pc    

         splice_Count    

         commonSampleIds = tissuewise_extractGenotype.out.genoVcfData_ch


 }