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
          tuple val(sample_id), file ("${sample_id}.junc"), emit: junc_ch 
        

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
         tuple val(sample_id), file ("${sample_id}.junc"), emit: junc_ch
        

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
          //file ("${sample_id}.junc") 

          tuple val(sample_id), file ("${sample_id}.junc"), emit: junc_ch
        

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
     
   tag "junc_input on $sample_id"
   input:
      

    tuple val(sample_id),file (juncfile)

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
        //file (leafcutter_cluster) 
        //file (phenotype_table)
        val(intron_min_len)
        val(inton_max_len)
        val(pheno_pcs) 
        

        output:
        
        path("output_perind.counts.gz.PCs"), emit: splicepcs

        path("Splice_count_Matrices_filtered.tsv"), emit: spliceCounts


        script:

        // -m 50 -l 500000 adopted from script 5 cGTEx paper DOI: 10.1038/s41588-022-01153-5

        //NOTE: User can change the deault number of PCs in the command below python ${phenotype_table} -p X

        /* leafcutter_cluster_Bovine_regtools.py will cluster together the introns fond in the junc files listed in junction_files.txt (created below in workflow), 
        requiring 50 split reads supporting each cluster and allowing introns of up to 500kb */

        /* prepare_phenotype_table.py script a) calculates intron excision ratios 
          b) filter out introns used in less than 40% of individuals or with almost no variation 
          c) output these ratios as gzipped txt files along with a user-specified number of PCs*/

        """   

        python $projectDir/bin/leafcutter/scripts/leafcutter_cluster_Bovine_regtools.py -j ${junc_file} -m $intron_min_len -l $inton_max_len  -o output
        
        python $projectDir/bin/leafcutter/scripts/prepare_phenotype_table.py -p $pheno_pcs output_perind.counts.gz

        awk 'FNR>1 || NR==1' output_perind.counts.gz.qqnorm_chr* > mergedCounts_tmp.txt

        cat mergedCounts_tmp.txt | awk -F'\t' 'NR>1{print \$1,\$2,\$3,\$4}' | awk -F'_' '{print \$1"_"\$2"\t"\$3}' | awk -v OFS='\t' '{print \$1,\$2,\$3,\$1":"\$2"-"\$3,\$4,\$5}' | sed '1i#Chr\tstart\tend\tpid\tgid\tstrand' > header.txt

        cat mergedCounts_tmp.txt | awk  '{for(i=5;i<=NF;i++) printf \$i"\t"; print ""}'  > counts.txt

        paste -d '\t' header.txt counts.txt | awk -v OFS='\t' 'NR==1 {print}' > xx

        paste -d'\t'  header.txt counts.txt |  awk -v OFS='\t' 'NR>1 {print}' | sort -k1,1d -k2,2n -k3,3n > xy

        cat xx xy | sed 's/#Chr/Chr/g' > Splice_count_Matrices_filtered.tsv


        """
}




process filter_Splicecounts {
        tag "on chromosome ${chr}"
        publishDir "${params.outdir}/Splicecounts_filtered", mode: 'copy'

        input:
        file (splcieCounts)
        
        output:
        

        tuple val(chr), file("fil_splcieCounts${chr}"), emit: spliceCounts_filt


        script:

        //prefix = splcieCounts[0].toString() - ~/(_R1)?(\.clean)?(\.fastq)?(\.gz)?$/

        """   

        cp $splcieCounts fil_splcieCounts${chr}

        """
}