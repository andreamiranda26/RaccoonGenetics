#!/bin/bash

#SBATCH --job-name=Samtools-idx
#SBATCH --output=logs/Samtools-idx-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=logs/Samtools-idx-%j.err           # Error file. %j is replaced with job ID
#SBATCH --cpus-per-task=20
#SBATCH --time=5-00
#SBATCH --mem=80G
#SBATCH -p jrw0107_std,general,nova

source $HOME/miniforge3/bin/activate Samtools

ROOT_DIR="/home/azm0272/Raccoons/Raccoon_analysis"

# Define variables
GENOME="GCA_028646535.1_pl-1k_genomic.fna"
GENOME_FILE="${ROOT_DIR}/data/reference/${GENOME}"
BAM="${ROOT_DIR}/data/bwamem2/Raccoon-Merged.bam"
BAI="${BAM}.bai"

samtools index -@ 20 --bai --output $BAI $BAM
