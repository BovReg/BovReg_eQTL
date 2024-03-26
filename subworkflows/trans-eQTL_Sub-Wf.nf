/*
========================================================================================
    Include Modules
========================================================================================
*/


include { transcriptCountsPCA_trans } from '../modules_dsl2/PhenotypePCs_trans-eQTL'

include { geneCountsPCA_trans} from '../modules_dsl2/PhenotypePCs_trans-eQTL'

include { rnaSplicePCS_trans } from '../modules_dsl2/PhenotypePCs_trans-eQTL'

include { trans_eQTL_nominal} from '../modules_dsl2/trans-eQTL'

include { trans_eQTL_permu} from '../modules_dsl2/trans-eQTL'

include { trans_eQTL_FDR} from '../modules_dsl2/trans-eQTL'



workflow transeQTL_workflow {

   take:

       genotypeStratificationPCA

       genotype_input_ch

       phenotypeTranscript_ch

       phenotypeGene_ch

       splicepcs_ch

       phenotypeSplice_ch

       ch_FDR_trans

       phenotype_PCs_trans

       threshold_trans_ch

       permutations_trans_ch

   main:

   /* phenotypes for trans-eQTL: Extracts the phenotypes from all chromosomes to perform trans-eQTL for each variant */

   transcriptCountsPCA_trans(phenotypeTranscript_ch.collect(),genotype_input_ch.join(genotypeStratificationPCA),phenotype_PCs_trans)

   geneCountsPCA_trans(phenotypeGene_ch.collect(), genotype_input_ch.join(genotypeStratificationPCA),phenotype_PCs_trans)

   rnaSplicePCS_trans(splicepcs_ch.collect(),phenotypeSplice_ch.collect(), genotype_input_ch.join(genotypeStratificationPCA))

   /* trans-eQTL nominal, permuataion and FDR */

   qtlmap_trans_geneCount = geneCountsPCA_trans.out.phenotype_gene_Bed_ch.join(geneCountsPCA_trans.out.cis_covariates_gene_ch).join(geneCountsPCA_trans.out.genotypeQTL_gene_ch)

   qtlmap_trans_transcriptCount = transcriptCountsPCA_trans.out.phenotype_transcript_Bed_ch.join(transcriptCountsPCA_trans.out.cis_covariates_transcript_ch).join(transcriptCountsPCA_trans.out.genotypeQTL_transcript_ch)
   
   qtlmap_trans_spliceCount = rnaSplicePCS_trans.out.phenotype_splice_Bed_ch.join(rnaSplicePCS_trans.out.cis_covariates_splice_ch).join(rnaSplicePCS_trans.out.genotypeQTL_splice_ch)

   trans_eQTL_nominal(qtlmap_trans_geneCount.join(qtlmap_trans_transcriptCount).join(qtlmap_trans_spliceCount), threshold_trans_ch)

   trans_eQTL_permu(qtlmap_trans_geneCount.join(qtlmap_trans_transcriptCount).join(qtlmap_trans_spliceCount), permutations_trans_ch,threshold_trans_ch )
   
   trans_eQTL_FDR(trans_eQTL_nominal.out.gene_trans_nominal_hits_ch.join(trans_eQTL_permu.out.gene_trans_permu_hits_ch).join(trans_eQTL_nominal.out.transcript_trans_nominal_hits_ch).join(trans_eQTL_permu.out.transcript_trans_permu_hits_ch).join(trans_eQTL_nominal.out.splice_trans_nominal_hits_ch).join(trans_eQTL_permu.out.splice_trans_permu_hits_ch),ch_FDR_trans)

  }




