/* This workflow takes the input RNAseq reads which are single-end and having stranded library fr-firststrand  */


/*
========================================================================================
    Include Modules
========================================================================================
*/


include { fastqc_SE } from '../modules_dsl2/QualityCheckandQualityControl'

include { trimmomatic_singleEnd} from '../modules_dsl2/QualityCheckandQualityControl'

include { starAlign_singleEnd_GeneTranscript } from '../modules_dsl2/Alignment'

include { starAlign_singleEnd_Splicing } from '../modules_dsl2/Alignment'

include {samtools_singleEnd_forleafcutter } from '../modules_dsl2/Indexing'



workflow SINGLE_END_READS{

        take:

          singelend_reads

          trimmjar_file

          adapter_file

          gtf_file

          star_index

          anchor_length

          intron_length_minimum

          intron_length_maximum

          stringtie_ip_di_ch

          leafcutter_ip_di_ch

          leafcutter_ip_idx_ch

        main:



         if ( params.bamFiles_input ) {
            
            stringtie_ip = stringtie_ip_di_ch

            leafcutter_bamtojunc_ip = leafcutter_ip_di_ch.join(leafcutter_ip_idx_ch)
           
         } 

         else{
         fastqc_SE(singelend_reads)

         trimmomatic_singleEnd(singelend_reads,trimmjar_file,adapter_file)

         starAlign_singleEnd_GeneTranscript(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         starAlign_singleEnd_Splicing(trimmomatic_singleEnd.out.forward_strand_trim ,gtf_ARS.collect(),star_index.collect())

         samtools_singleEnd_forleafcutter(starAlign_singleEnd_Splicing.out.samindex_ch)
    
         stringtieQuant_library_firststranded(starAlign_singleEnd_GeneTranscript.out.stringtie_ch,gtf_ARS.collect())

         leafcutter_bamtojunc_library_firststranded(samtools_singleEnd_forleafcutter.out.tosingelEndleafcutter_ch,anchor_length,intron_length_minimum,intron_length_maximum) 
         

         }
       
         emit:
    
         stringtie_ip

         leafcutter_bamtojunc_ip

}