

 /*
                                  STEP 4 - Read indexing with samtools and the output BAM files are supplied
                                            to leafcutter for extractiong splice junctions
                                            */

process samtools_merge_forStringtie {
        tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort", mode: 'copy'

        input:
         //tuple val(sample_id), file(samindex)
        tuple val(sample_id), file(paired),file(unpairedR1),file(unpairedR2)

        output:

         tuple val(sample_id), file("${sample_id}_merged_sorted.bam"), emit:tostringtie_ch

        script:
        
 

        """ 
        samtools merge ${sample_id}_merged.bam ${paired} ${unpairedR1} ${unpairedR2}
        
        samtools sort ${sample_id}_merged.bam -o ${sample_id}_merged_sorted.bam --threads 18
        
        rm ${sample_id}_merged.bam

        """
   
}




process samtools_merge_forleafcutter  {
        tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort", mode: 'copy'

        input:

        tuple val(sample_id), file(pairedLeafcutter),file(unpairedLeafR1),file(unpairedLeafR2)

        output:

        tuple val(sample_id), file("${sample_id}_leafcutter_merged_sorted.bam"),file("${sample_id}_leafcutter_merged_sorted.bai"), emit:toleafcutter_ch


        script:
        

          
        """ 
        samtools merge ${sample_id}_merged.bam ${pairedLeafcutter} ${unpairedLeafR1} ${unpairedLeafR2}

        samtools sort ${sample_id}_merged.bam -o ${sample_id}_leafcutter_merged_sorted --threads 18
        samtools index  ${sample_id}_leafcutter_merged_sorted -@ 18
        mv ${sample_id}_leafcutter_merged_sorted ${sample_id}_leafcutter_merged_sorted.bam

        rm ${sample_id}_merged.bam

        """

}




process samtools_singleEnd_forleafcutter  {
        tag "samtools index on $sample_id"
        container 'biocontainers/samtools:v1.9-4-deb_cv1'
        publishDir "${params.outdir}/samSort", mode: 'copy'

        input:
        tuple val(sample_id), file(samindex)
        

        output:

        tuple val(sample_id), file("${sample_id}_leafcutter_sorted.bam"),file("${sample_id}_leafcutter_sorted.bai"), emit:tosingelEndleafcutter_ch


        script:
        

         
        """ 
        samtools sort ${samindex} -o ${sample_id}_leafcutter_sorted --threads 18
        samtools index  ${sample_id}_leafcutter_sorted -@ 18
        mv ${sample_id}_leafcutter_sorted ${sample_id}_leafcutter_sorted.bam

        """

}


