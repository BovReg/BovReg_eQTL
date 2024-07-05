## Demo data

- The genotype should be provided in [vcf](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format, a demo data for test run can be downloaded here: Demo_genotype_data [download](https://zenodo.org/records/10997393/files/Demo_genotype_BovReg.tar.gz?download=1) and a [Demodata/Geno_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Geno_input.tsv) file was provided with the file paths of different chromosomes present in the demo data. ( Note: For the trial run the downloaded folder containing all the vcf.gz files "Demo_genotype_BovReg.tar.gz" should be saved and uncompressed in the folder [Demo_genotype_BovReg](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Demo_genotype_BovReg)).

- The Phenotype data can be provided in any of the following formats and all the input files paths should be provides in a tsv file. An example tsv files with demodata is available in the folder [Demodata](https://github.com/BovReg/BovReg_eQTL/tree/main/Demodata).

   1. Raw data (RNAseq expression data in fastq format): Demo_data [download](https://zenodo.org/records/7949616/files/Demo_RNAseqData_BovReg.tar.gz?download=1) and a [Demodata/fastq_paired_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/fasta_paired_input.tsv) file was provided with the file paths for different fastq samples present in the demo data. ( Note: For the trial run the downloaded folder containing the corresponding fastq files "Demo_RNAseqData_BovReg.tar.gz" should be saved and uncompressed in the folder [Demo_RNAseqFastq](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Demo_RNAseqFastq)) . 

   2. Aligned reads (RNAseq expression data in bam format): Demo_data [download](https://zenodo.org/records/7950181/files/Demo_RNAseqBam_BovReg.tar.gz?download=1) and a [Demodata/Bam_input.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Bam_input.tsv) file was provided with the file paths for different bam samples present in the demo data. ( Note: For the trial run the downloaded folder containing the corresponding bam files "Demo_RNAseqBam_BovReg.tar.gz" should be saved and uncompressed in the folder [Demo_RNAseqBam](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Demo_RNAseqBam)) .

   4. Expression counts across samples for gene level, transcript level and splicing counts (expression count matrices as text file): Demo_data [Demodata/Count_matrices.tsv](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/Count_matrices.tsv). 
   

- The genotype-phenotype corresponding samples information should be provided as text file: [Demo_data/RNA_WGS_CorresID_BovReg.txt](https://github.com/BovReg/BovReg_eQTL/blob/main/Demodata/RNA_WGS_CorresID_BovReg.txt).




