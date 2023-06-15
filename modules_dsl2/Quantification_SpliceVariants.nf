                            /*
                             * STEP 6 -Quantification of RNA splicing Leafcutter (intron excisions)
                             */

/*
 * Leafcutter bam to junctions
 * The junction lines with "?" were removed as they are generating error in the cluster perperation
 */
 
   //leafcutter_junc_ip = leafcutter_junc_bam_ch.merge(leafcutter_junc_bai_ch)
process leafcutter_bamtojunc_library_firststranded {
        tag "leafcutter_junc on $sample_id"
        container 'griffithlab/regtools:release-0.6.0'
        publishDir "${params.outdir}/leafcutter_bamtojunc", mode: 'copy'

        input:
        
        tuple val(sample_id), file (samtool_bam), file(samtool_bai)
        val(anchor_len)
        val(intron_min_len)
        val(inton_max_len)

         
        
        output:
          file ("${sample_id}.junc") 
        

        script:


        """
        regtools junctions extract -s 1 -a $anchor_len -m $intron_min_len -M $inton_max_len ${samtool_bam} -o ${sample_id}.junc

        grep -v "?" ${sample_id}.junc > ${sample_id}_Mod.junc

        mv ${sample_id}_Mod.junc ${sample_id}.junc

        """
}



process leafcutter_bamtojunc_library_secondstranded {
        tag "leafcutter_junc on $sample_id"
        container 'griffithlab/regtools:release-0.6.0'
        publishDir "${params.outdir}/leafcutter_bamtojunc", mode: 'copy'

        input:
        
        tuple val(sample_id), file (samtool_bam), file(samtool_bai)
        val(anchor_len)
        val(intron_min_len)
        val(inton_max_len)         
        
        output:
          file ("${sample_id}.junc") 
        

        script:


        """
        regtools junctions extract -s 2 -a $anchor_len -m $intron_min_len -M $inton_max_len ${samtool_bam} -o ${sample_id}.junc

        grep -v "?" ${sample_id}.junc > ${sample_id}_Mod.junc

        mv ${sample_id}_Mod.junc ${sample_id}.junc

        """
}



process leafcutter_bamtojunc_library_unstranded {
        tag "leafcutter_junc on $sample_id"
        container 'griffithlab/regtools:release-0.6.0'
        publishDir "${params.outdir}/leafcutter_bamtojunc", mode: 'copy'

        input:
        
        tuple val(sample_id), file (samtool_bam), file(samtool_bai)
        val(anchor_len)
        val(intron_min_len)
        val(inton_max_len)         
        
        output:
          file ("${sample_id}.junc") 
        

        script:


        """
        regtools junctions extract -s 0 -a $anchor_len -m $intron_min_len -M $inton_max_len ${samtool_bam} -o ${sample_id}.junc

        grep -v "?" ${sample_id}.junc > ${sample_id}_Mod.junc

        mv ${sample_id}_Mod.junc ${sample_id}.junc

        """
}



/* Extract jun files from different paths and for making the nexflow to recorgnise the path in 'work' folder

*/


process rnaspliceFilterJunc {
     

   input:
      

    file (juncfile)

   output:

     file ("comm_geno_${juncfile}")

  
   script:

    juncfilename = juncfile[0].toString() 
    
  // grep -q ${juncfilename} samplelist && echo "\$(cat juncfile)"  > ${juncfilename}.junc
     // for f in `cat samplelistjunc`; do cp \${f}  comm_\${f}; done;

    """  
         cp ${juncfile} comm_geno_${juncfile}
         
    """
}


/*
 * Leafcutter quantification step
   The leafcutter_cluster file by default clusters for 22 chromsomes (human genome) and for the cattle genome a minor tweak ie required
    The range for chromsome was changes from 1-23 to 1-30 in the new script "$baseDir/leafcutter/scripts/leafcutter_cluster_Bovine_regtools.py"
    This was used for this process to cover all the bovine autosomal chromosomes.

 */

process leafcutter_cluster_junctions {
        tag "leafcutter_cluster_Junctions"
        container 'faizanbashir/python-datascience:2.7'
        publishDir "${params.outdir}/leafcutter_cluster", mode: 'copy'

        input:
        file (junc_file) 
        file (leafcutter_cluster) 
        file (phenotype_table)
        val(intron_min_len)
        val(inton_max_len)
        val(pheno_pcs) 
        output:
        
        file ("*") 

        script:

        // -m 50 -l 500000 adopted from script 5 cGTEx paper DOI: 10.1038/s41588-022-01153-5

        //NOTE: User can change the deault number of PCs in the command below python ${phenotype_table} -p X
        """   
   
        python ${leafcutter_cluster} -j ${junc_file} -m $intron_min_len -l $inton_max_len  -o output
        
        python ${phenotype_table} -p $pheno_pcs output_perind.counts.gz
        """
}

