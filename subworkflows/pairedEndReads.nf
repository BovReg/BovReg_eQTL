/* This workflow takes the input RNAseq reads which are paired-end and having stranded library fr-firststrand  */


/*
========================================================================================
    Include Modules
========================================================================================
*/

include { fastqc_PE } from '../modules_dsl2/QualityCheckandQualityControl'

include { trimmomatic } from '../modules_dsl2/QualityCheckandQualityControl'

include { starAlign_GeneTranscript } from '../modules_dsl2/Alignment'

include { starAlign_Splicing } from '../modules_dsl2/Alignment'

include { starAlign_unpaired } from '../modules_dsl2/Alignment'

include {samtools_merge_forStringtie} from '../modules_dsl2/Indexing'

include {samtools_merge_forleafcutter } from '../modules_dsl2/Indexing'




params.bamFiles_input  = false



workflow PAIREDEND_END_READS {

        take:

          paired_reads

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
             fastqc_PE(paired_reads)

             trimmomatic(paired_reads,trimmjar_file,adapter_file)

             starAlign_GeneTranscript(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

             starAlign_Splicing(trimmomatic.out.forward_strand_trim.join(trimmomatic.out.revese_strand_trim),gtf_file.collect(),star_index.collect())

             starAlign_unpaired(trimmomatic.out.forward_unpaired_trim.join(trimmomatic.out.revese_unpaired_trim),gtf_file.collect(),star_index.collect())
       
             samtools_merge_forStringtie(starAlign_GeneTranscript.out.tomerge_paired_ch.join(starAlign_unpaired.out.tomerge_unpairedR1_ch).join(starAlign_unpaired.out.tomerge_unpairedR2_ch))

             samtools_merge_forleafcutter(starAlign_Splicing.out.tomerge_paired_Leaf_ch.join(starAlign_unpaired.out.tomerge_unpairedLeafR1_ch).join(starAlign_unpaired.out.tomerge_unpairedLeafR2_ch))
            
             stringtie_ip = samtools_merge_forStringtie.out.tostringtie_ch

             leafcutter_bamtojunc_ip = samtools_merge_forleafcutter.out.toleafcutter_ch

            }


         emit:
    
         stringtie_ip

         leafcutter_bamtojunc_ip


}
