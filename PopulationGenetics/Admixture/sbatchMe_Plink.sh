#!/bin/bash

#SBATCH --job-name=Plink
#SBATCH --output=logs/Plink-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=logs/Plink-%j.err           # Error file. %j is replaced with job ID
#SBATCH --cpus-per-task=2
#SBATCH --time=5-00
#SBATCH --mem=10G
#SBATCH -p jrw0107_std,general,nova_20

source $HOME/miniforge3/bin/activate Admixture

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial" 

# Define variables
IN_VCF="${ROOT_DIR}/data/vcf/Raccoon.filtered.vcf.gz"
OUTPUT="${ROOT_DIR}/analysis/admixture/admixture"

# Generate the input file in plink format
plink --vcf $IN_VCF --make-bed --out $OUTPUT --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $OUTPUT.bim > $OUTPUT.bim.tmp
mv $OUTPUT.bim.tmp $OUTPUT.bim