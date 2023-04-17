

/*
process TranscriptGTFfiles:
       -This process takes the Stringtie gtf files as input_files
       -It extract transcripts from each gtf and process the headers
*/

process transcriptGTFfiles {
  tag "$gtf_input"
  input:
    
   // tuple val(sample_id), file (gene_input)
    
    file (gtf_input)  

  output:

   file ("${sample_id}_transcriptCounts.txt") 

  script:

  sample_id = gtf_input[0].toString() - ~/(\_stringtie.gtf)?$/
  """
  awk -v OFS='\t' 'NR==1 {gsub(/_stringtie.gtf/,"",\$11) ;print \$1="Chr",\$2="start",\$3="end",\$4="pid",\$5="gid",\$6="strand",\$7="temp" }; \
  NR>1{if(\$3=="transcript" && \$19=="TPM") print \$1,\$4,\$5,\$12,\$10,\$7,\$20 ; \
  else if(\$3=="transcript" && \$17=="TPM") print \$1,\$4,\$5,\$12,\$10,\$7,\$18}' ${gtf_input} \
   | sed 's/"//g;s/;//g' | sed 's/temp/${sample_id}/g' > "${sample_id}_transcriptCounts.txt"


  """
}




/*
 process MergeNormalizedTranscriptCountMatrices:
    - The TPM normalized transcript counts for each sample were merged into a single file to perform PCA for eQTL detection.
    -The input file comes form Stringtie .gtf files
    -The TPM normalization was already performed by Stringtie
    -The positions, with atleast ten percent samples have TPM values greater than 1 will be filtered out

*/

process mergeNormalizedTranscriptCountMatrices {
   tag "$input_files"
   publishDir "${params.outdir}/Count_Matrices", mode:'copy'
   container 'quay.io/biocontainers/csvtk:0.21.0--0'

   input:
    file (input_files) 

  output:
    file ("Transcript_count_Matrices_filtered.tsv")
  script:

  def single = input_files instanceof Path ? 1 : input_files.size()
  def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "pid,gid,start,end,Chr,strand"'
  """
  $merge $input_files >  Transcript_count_Matrices.tsv
  awk 'NR==1;NR>1   {z=0; for (i=1; i<=NF; ++i) if (\$i < 1 && ++z> 0.9*(NF-7) ) next; print }'  Transcript_count_Matrices.tsv \
   > Transcript_count_Matrices_filtered.tsv
  """

}
