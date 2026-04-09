#!/bin/bash
#SBATCH --job-name=join_landcover
#SBATCH --output=logs/join_landcover-%j.out
#SBATCH --error=logs/join_landcover-%j.err
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

# activate conda environment
source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
SCRIPT="${ROOT_DIR}/scripts/PopulationGenetics/join_labels_landcover.R"

mkdir -p ${ROOT_DIR}/scripts/PopulationGenetics/logs

echo "Starting label + landcover join"
date

Rscript ${SCRIPT}

echo "Finished"
date