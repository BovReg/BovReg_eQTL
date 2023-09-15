

# Built image from the commands from "Dockerfile"

docker build -t praveen/rpackages-eqtl .


# Run docker r studio

docker run -it --rm --user rstudio praveen/rpackages-eqtl R

# Run bash script

docker run -it --rm --user rstudio praveen/rpackages-eqtl bash

# Mount directory and set working directory
docker run -it --rm -v $(pwd):/disk2/chitneedi/NextFlow/eQTL_Genotype_vcf/liver_Test_Genotypes -w /disk2/chitneedi/NextFlow/eQTL_Genotype_vcf/liver_Test_Genotypes praveen/rpackages-eqtl R

