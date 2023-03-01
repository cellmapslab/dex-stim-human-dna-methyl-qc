#!/bin/bash
#SBATCH --job-name="getBMIQbetas"
#SBATCH --part=pe
#SBATCH --mem=200GB
#SBATCH --output=getBMIQbetas.out

module load R

SRC_DIR=~mpip/code/dex-stim-human-dna-methyl-qc/

Rscript --vanilla $SRC_DIR/03_normalization/05_getBMIQbetas.R $SRC_DIR/input_parameters.csv