# BovReg_eQTL
Repository
This is the BovReg_eQTL analysis.
- This repository provide the nextflow scripts and demo data to test eQTL analysis and run with large datasets.
- User should have a linux environment with the lastest nextflow and docker installed to run these scripts.
- Care should be taken with the input and output paths while running with new data.
- The numbering of the scripts is the order in which the scripts are executed with nextflow.
- To run a test analysis, demodata is available in the folder with the same name.
- The script 1,3 and 4 contain json files for calling parmaeters, whcih can run using the command 
        nextflow run script.nf -params-file file.json 
