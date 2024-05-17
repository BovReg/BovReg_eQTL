
/*
========================================================================================
    Include Modules
========================================================================================
*/


include { transcriptCountsPCA } from '../modules_dsl2/PhenotypePCs'

include { geneCountsPCA } from '../modules_dsl2/PhenotypePCs'

include { rnaSplicePCS } from '../modules_dsl2/PhenotypePCs'

include { cisQTL_nominal } from '../modules_dsl2/cis-eQTL'

include { cisQTL_permutation } from '../modules_dsl2/cis-eQTL'

include { cisQTL_conditional} from '../modules_dsl2/cis-eQTL'


workflow ciseQTL_workflow {

     take:

       genotypeStratPCA

       genotype_input_ch

       phenotypeTranscript_ch

       phenotypeGene_ch

       splicepcs_ch

       phenotypeSplice_ch

       phenotype_PCs_cis
   
       nominal_cis_ch

       permutations_cis_ch

       fdr_rate_cis_ch

    main:
  
   /* phenotypes for cis-eQTL: Extracts the phenotypes from each chromosomes to perform chromosome wise cis-eQTL */

   transcriptCountsPCA(phenotypeTranscript_ch.collect(),genotype_input_ch.join(genotypeStratPCA),phenotype_PCs_cis)

   geneCountsPCA(phenotypeGene_ch.collect(), genotype_input_ch.join(genotypeStratPCA),phenotype_PCs_cis)
   
   rnaSplicePCS(phenotypeSplice_ch.collect(), genotype_input_ch.join(genotypeStratPCA),splicepcs_ch.collect())

  /*cis-eQTL nominal and permuataion*/

   qtlmap_cis_geneCount = geneCountsPCA.out.phenotype_gene_Bed_ch.join(geneCountsPCA.out.cis_covariates_gene_ch).join(geneCountsPCA.out.genotypeQTL_gene_ch)

   qtlmap_cis_transcriptCount = transcriptCountsPCA.out.phenotype_transcript_Bed_ch.join(transcriptCountsPCA.out.cis_covariates_transcript_ch).join(transcriptCountsPCA.out.genotypeQTL_transcript_ch)
   
   qtlmap_cis_spliceCount = rnaSplicePCS.out.phenotype_splice_Bed_ch.join(rnaSplicePCS.out.cis_covariates_splice_ch).join(rnaSplicePCS.out.genotypeQTL_splice_ch)

   cisQTL_nominal(qtlmap_cis_geneCount.join(qtlmap_cis_transcriptCount).join(qtlmap_cis_spliceCount),nominal_cis_ch)
   
   cisQTL_permutation(qtlmap_cis_geneCount.join(qtlmap_cis_transcriptCount).join(qtlmap_cis_spliceCount),permutations_cis_ch)

   /*cis-eQTL conditional */
/*
   qtlmap_cis_geneCount_cond = qtlmap_cis_geneCount.join(cisQTL_permutation.out.cis_permu_geneResults_ch)

   qtlmap_cis_transcriptCount_cond = qtlmap_cis_transcriptCount.join(cisQTL_permutation.out.cis_permu_transcritResults_ch)

   qtlmap_cis_spliceCount_cond = qtlmap_cis_spliceCount.join(cisQTL_permutation.out.cis_permu_splicingResults_ch)

   cisQTL_conditional(qtlmap_cis_geneCount_cond.join(qtlmap_cis_transcriptCount_cond).join(qtlmap_cis_spliceCount_cond),fdr_rate_cis_ch)
*/
}