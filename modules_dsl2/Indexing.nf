

/* Indexing of bam files was not required for stringtie quntification */


process samtools_merge_forStringtie{
        //tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort", mode: 'copy'

        input:
         //tuple val(sample_id), file(samindex)
        tuple val(sample_id), file(paired),file(unpairedR1),file(unpairedR2)

        output:
         //tuple  val(sample_id), file("*.bam"),file("*.bai") 
         tuple val(sample_id), file("${samplename}_merged_sorted.bam")

        script:
        
        samplename = paired[0].toString() - ~/(_paired)?(_Aligned)?(\.sortedByCoord)?(\.out)?(\.bam)?$/

        // after alignment the paired and unpaired reads were merged and sorted for stringtie input

        """ 
        samtools merge ${samplename}_merged.bam ${paired} ${unpairedR1} ${unpairedR2}
        
        samtools sort ${samplename}_merged.bam -o ${samplename}_merged_sorted.bam --threads 18
        

        """
   //        samtools view -h -F 4 ${leafcutter_bam}.sorted  | samtools view -bS > ${leafcutter_bam}.tmp
}






/* Indexing of bam files was not required for stringtie quntification */


process samtools_merge_forStringtie_with_rf {
        //tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort_rf", mode: 'copy'

        input:
         //tuple val(sample_id), file(samindex)
        tuple val(sample_id), file(paired),file(unpairedR1),file(unpairedR2)

        output:
         //tuple  val(sample_id), file("*.bam"),file("*.bai") 
         tuple val(sample_id), file("${sample_id}_merged_sorted.bam"), emit:tostringtie_with_rf

        script:
        
        //samplename = paired[0].toString() - ~/(_paired)?(_Aligned)?(\.sortedByCoord)?(\.out)?(\.bam)?$/

        // after alignment the paired and unpaired reads were merged and sorted for stringtie input

        """ 
        samtools merge ${sample_id}_merged.bam ${paired} ${unpairedR1} ${unpairedR2}
        
        samtools sort ${sample_id}_merged.bam -o ${sample_id}_merged_sorted.bam --threads 18
        

        """
   //        samtools view -h -F 4 ${leafcutter_bam}.sorted  | samtools view -bS > ${leafcutter_bam}.tmp
}






process samtools_merge_forStringtie_without_rf {
        //tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort_norf", mode: 'copy'

        input:
         
        tuple val(sample_id), file(paired),file(unpairedR1),file(unpairedR2)

        output:
         tuple val(sample_id), file("${sample_id}_merged_sorted.bam"), emit:tostringtie_without_rf

        script:
        
        //samplename = paired[0].toString() - ~/(_paired)?(_Aligned)?(\.sortedByCoord)?(\.out)?(\.bam)?$/

        // after alignment the paired and unpaired reads were merged and sorted for stringtie input

        """ 
        samtools merge ${sample_id}_merged.bam ${paired} ${unpairedR1} ${unpairedR2}
        
        samtools sort ${sample_id}_merged.bam -o ${sample_id}_merged_sorted.bam --threads 18
        

        """
   //        samtools view -h -F 4 ${leafcutter_bam}.sorted  | samtools view -bS > ${leafcutter_bam}.tmp
}



 /*
                                  STEP 4 - Read indexing with samtools and the output BAM files are supplied
                                            to leafcutter for extractiong splice junctions
                                            */

process samtools_merge_forleafcutter  {
        //tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort", mode: 'copy'

        input:
         //tuple val(sample_id), file(samindex)
        tuple val(sample_id), file(pairedLeafcutter),file(unpairedLeafR1),file(unpairedLeafR2)

        output:
         //tuple  val(sample_id), file("*.bam"),file("*.bai") 
        tuple file("${samplename}_leafcutter_merged_sorted.bam"),file("${samplename}_leafcutter_merged_sorted.bai")

        //samtools view $leafcutter_bam | python $filter_cs | $sam2bed --use-RNA-strand - ${leafcutter_bam.baseName}.bed
        //$bed2junc ${leafcutter_bam.baseName}.bed ${leafcutter_bam.baseName}.junc

        script:
        
        samplename = pairedLeafcutter[0].toString() - ~/(_paired)?(_leafcutter)?(_Aligned)?(.out)?(.bam)?$/
        //  cp $leafcutter_bam ${leafcutter_bam}.tmp
        // Removing the reads unmapped to a chromosome is required for process leafcutter_cluster_junctions
        // The reads need to be sorted for indexing using samtools and it was not mentioned in leafcutter documentation (webpage)
        // after alignment the paired and unpaired reads were merged and sorted 
        // samtools sort -n has been used to sort the reads by name instead. Remove -n to sort by position, which is what is needed to prepare a BAM file for indexing with samtools index
          
        """ 
        samtools merge -f ${samplename}_merged.bam ${pairedLeafcutter} ${unpairedLeafR1} ${unpairedLeafR2}

        samtools sort ${samplename}_merged.bam -o ${samplename}_leafcutter_merged_sorted --threads 18
        samtools index  ${samplename}_leafcutter_merged_sorted -@ 18
        mv ${samplename}_leafcutter_merged_sorted ${samplename}_leafcutter_merged_sorted.bam

        """
   //        samtools view -h -F 4 ${leafcutter_bam}.sorted  | samtools view -bS > ${leafcutter_bam}.tmp
}







process samtools_singleEnd_forleafcutter  {
        //tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort", mode: 'copy'

        input:
        tuple val(sample_id), file(samindex)
        

        output:
        //tuple  val(sample_id), file("*.bam"),file("*.bai") 
        tuple file("${samplename}_leafcutter_sorted.bam"),file("${samplename}_leafcutter_sorted.bai")

        //samtools view $leafcutter_bam | python $filter_cs | $sam2bed --use-RNA-strand - ${leafcutter_bam.baseName}.bed
        //$bed2junc ${leafcutter_bam.baseName}.bed ${leafcutter_bam.baseName}.junc

        script:
        
        samplename = samindex[0].toString() - ~/(_leafcutter)?(_Aligned)?(.out)?(.bam)?$/
        //  cp $leafcutter_bam ${leafcutter_bam}.tmp
        // Removing the reads unmapped to a chromosome is required for process leafcutter_cluster_junctions
        // The reads need to be sorted for indexing using samtools and it was not mentioned in leafcutter documentation (webpage)
        // samtools sort -n has been used to sort the reads by name instead. Remove -n to sort by position, which is what is needed to prepare a BAM file for indexing with samtools index

        // samtools sort ${samindex} -o ${sample_id}_sorted --threads 18
        // samtools index  ${sample_id}_sorted -@ 18
        // mv ${sample_id}_sorted ${sample_id}_sorted.bam
         
        """ 
        samtools sort ${samindex} -o ${samplename}_leafcutter_sorted --threads 18
        samtools index  ${samplename}_leafcutter_sorted -@ 18
        mv ${samplename}_leafcutter_sorted ${samplename}_leafcutter_sorted.bam

        """
   //        samtools view -h -F 4 ${leafcutter_bam}.sorted  | samtools view -bS > ${leafcutter_bam}.tmp
}
