
/*
process Genetsvfiles:
       -This process takes the Stringtie tsv files as input_files
       -It extract gene normalized counts from each tsv and process the headers
*/

process geneTSVfiles {
  tag " $gene_input"
  input:
  // tuple val(sample_id), file (gene_input) 
     file (gene_input)
  output:

   file ("${sample_id}_gene.txt") 

  script:

  sample_id = gene_input[0].toString() - ~/(\_stringtie.tsv)?$/
  """
  awk -v OFS='\t' 'NR==1 {print \$1="Chr",\$2="start",\$3="end",\$4="pid",\$5="gid",\$6="strand","${sample_id}" }; \
  NR>1{ print \$3,\$5,\$6,\$1,\$1,\$4,\$9 }' ${gene_input}  > ${sample_id}_gene.txt
  """
}

/*
 process MergeNormalizedTranscriptCountMatrices:
    - The TPM normalized transcript counts for each sample were merged into a single file to perform PCA for eQTL detection.
    -The input file comes form Stringtie .gtf files
    -The TPM normalization was already performed by Stringtie
    -The positions with atleast ten percent samples have TPM values greater than 1 will be filtered out

*/

process mergeNormalizedGeneCountMatrices {

   publishDir "${params.outdir}/Count_Matrices", mode:'copy'
   container 'quay.io/biocontainers/csvtk:0.21.0--0'

   input:
    file (input_files) 

  output:
  
  file ("Gene_count_Matrices_filtered.tsv")
  script:

  def single = input_files instanceof Path ? 1 : input_files.size()
  def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "pid,gid,start,end,Chr,strand"'
  """
  $merge $input_files >  Gene_count_Matrices.tsv
  awk 'NR==1;NR>1   {z=0; for (i=1; i<=NF; ++i) if (\$i < 1 && ++z> 0.9*(NF-7) ) next; print }' Gene_count_Matrices.tsv \
   > Gene_count_Matrices_filtered.tsv
  """

}