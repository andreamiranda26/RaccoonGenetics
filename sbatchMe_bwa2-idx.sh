#!/bin/bash
#to be able to map to the genome, you have to index first and then the runMe-bwa2_mem.sh runs the mapping 

#SBATCH --job-name=BWA2-Index
#SBATCH --output=logs/BWA2-Index-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=logs/BWA2-Index-%j.err           # Error file. %j is replaced with job ID
#SBATCH --cpus-per-task=5
#SBATCH --time=5-00
#SBATCH --mem=100G
#SBATCH -p jrw0107_std,general,nova

source $HOME/miniforge3/bin/activate Bwa # bwa-mem2 2.2.1

ROOT_DIR="/home/azm0272/Raccoons/Raccoon_analysis"

# Define variables
GENOME="GCA_028646535.1_pl-1k_genomic.fna"
GENOME_FILE="${ROOT_DIR}/data/reference/${GENOME}"

# Run BWA-MEM2
bwa-mem2 index ${GENOME_FILE}
