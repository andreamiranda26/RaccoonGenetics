#!/bin/bash
#SBATCH --job-name=NLCD_county
#SBATCH --output=logs/NLCD-%j.out
#SBATCH --error=logs/NLCD-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00
#SBATCH --mem=16G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
SCRIPT="${ROOT_DIR}/scripts/PopulationGenetics/extract_nlcd_by_county.R"

mkdir -p ${ROOT_DIR}/scripts/PopulationGenetics/logs

echo "Starting NLCD county extraction"
date

Rscript ${SCRIPT}

echo "Finished"
date