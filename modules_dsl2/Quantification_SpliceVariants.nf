                           /*
                             * STEP 6 -Quantification of RNA splicing Leafcutter (intron excisions)
                             */

/*
 * Leafcutter bam to junctions
 * The junction lines with "?" were removed as they are generating error in the cluster perperation
 */
 
   //leafcutter_junc_ip = leafcutter_junc_bam_ch.merge(leafcutter_junc_bai_ch)
process leafcutter_junctions {
        //tag "leafcutter_junc on $sample_id"
        container 'griffithlab/regtools:release-0.6.0'
        publishDir "${params.outdir}/leafcutter_bamtojunc", mode: 'copy'

        input:
        
        //tuple val(sample_id), file (samtool_bam), file(samtool_bim)
         
        tuple file (samtool_bam), file(samtool_bim)
        
        output:
         //file ("${sample_id}.junc") 

         file ("${samplename}.junc")
        

        script:
        samplename = samtool_bam[0].toString() - ~/(_leafcutter)?(_merged)?(_sorted)?(.bam)?$/

        
        //regtools junctions extract -s 2 -a 8 -m 50 -M 500000 ${samtool_bam} -o ${sample_id}.junc

        //grep -v "?" ${sample_id}.junc > ${sample_id}_Mod.junc

        //mv ${sample_id}_Mod.junc ${sample_id}.junc

        """
        regtools junctions extract -s 2 -a 8 -m 50 -M 500000 ${samtool_bam} -o ${samplename}.junc

        grep -v "?" ${samplename}.junc > ${samplename}_Mod.junc

        mv ${samplename}_Mod.junc ${samplename}.junc

        """
}


/* Extract jun files from different and check whether the correpoding genotypes are found
  
   if exits the junc file is ussed for clustering 

*/


process rnaspliceFilterJunc {
     container 'praveen/qtltools'

   input:
      

      file (juncfile)

   output:

     file ("comm_${juncfile}")

  
   script:

    juncfilename = juncfile[0].toString() 
    
  // grep -q ${juncfilename} samplelist && echo "\$(cat juncfile)"  > ${juncfilename}.junc
     // for f in `cat samplelistjunc`; do cp \${f}  comm_\${f}; done;

    """  
         cp ${juncfile} comm_${juncfile}
   
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
        output:
        
        file ("output_*") 

        script:

        // -m 50 -l 500000 adopted from script 5 cGTEx paper
        """   
   
        python ${leafcutter_cluster} -j ${juncfile} -m 50 -l 500000  -o output
        
        python ${phenotype_table} -p 15 output_perind.counts.gz
        """
}


process leafcutter_cluster_counts {
        tag "leafcutter_cluster_Junctions"
        container 'faizanbashir/python-datascience:2.7'
        publishDir "${params.outdir}/leafcutter_cluster", mode: 'copy'

        input:

        file (output_counts) 
        file (phenotype_table) 
        output:
        //file 'output*' into leafcutter_cluster_ch
        //file 'junction_files.tar.gz' into junction_files_ch
        file ("*") 

        script:

        // -m 50 -l 500000 adopted from script 5 cGTEx paper
        """
        
        python ${phenotype_table} -p 15 output_perind.counts.gz
        """
}