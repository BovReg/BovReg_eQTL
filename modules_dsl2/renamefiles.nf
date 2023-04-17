
/* The junc files from a tissue but from different sub-categopry with same sample ID were renamed  */

process newID_junc {

   tag "leafcutter_junc on $sample_id"
     
   input:
      
    tuple val(sampleID), file (juncfile_LL), file(juncfile_GM)

   output:

    path ("*"),  emit: junc_files

   script:

    juncfilename = juncfile[0].toString() 

    """  
         cp ${juncfile_LL} ${sampleID}_LL.junc

         
    """
}