#!/bin/bash

#SBATCH --job-name=Samtools-Merge
#SBATCH --output=logs/Samtools-Merge-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=logs/Samtools-Merge-%j.err           # Error file. %j is replaced with job ID
#SBATCH --cpus-per-task=30
#SBATCH --time=5-00
#SBATCH --mem=100G
#SBATCH -p jrw0107_std,general,nova

source $HOME/miniforge3/bin/activate Samtools

ROOT_DIR="/home/azm0272/Raccoons/Raccoon_analysis"

# Define variables
GENOME="GCA_028646535.1_pl-1k_genomic.fna"
GENOME_FILE="${ROOT_DIR}/data/reference/${GENOME}"
IN_BAMS="${ROOT_DIR}/data/bwamem2/*.bam"
OUT_BAM="${ROOT_DIR}/data/bwamem2/Raccoon-Merged.bam"

samtools merge --threads 30 --reference ${GENOME_FILE} -o ${OUT_BAM} ${IN_BAMS}
