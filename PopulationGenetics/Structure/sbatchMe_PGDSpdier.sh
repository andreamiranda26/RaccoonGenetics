#!/usr/bin/env bash

#SBATCH --job-name=PGDSpider
#SBATCH --cpus-per-task=1
#SBATCH --partition=jrw0107_std,general,nova_20
#SBATCH --time=1:00:00              # Run time (D-HH:MM:SS)
#SBATCH --output=logs/PGDSpider-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=logs/PGDSpider-%j.err           # Error file. %j is replaced with job ID

source $HOME/miniforge3/bin/activate Structure

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"

# Define variables

IN_VCF="${ROOT_DIR}/data/vcf/Raccoon.filtered.vcf"
OUTPUT="${ROOT_DIR}/analysis/structure/Raccoon.str"
SPID="${ROOT_DIR}/scripts/PopulationGenetics/Structure/VCF-STRUCTURE_dip.spid"

mkdir -p ${ROOT_DIR}/analysis/structure

PGDSpider2-cli -inputfile $IN_VCF -inputformat vcf -outputfile $OUTPUT -outputformat structure -spid $SPID