#!/bin/bash
#SBATCH --job-name=rabies_orv_map
#SBATCH --output=logs/rabies_orv_map-%j.out
#SBATCH --error=logs/rabies_orv_map-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

# ---- activate your conda env ----
source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
SCRIPT="${ROOT_DIR}/scripts/PopulationGenetics/plot_rabies_orv_choropleth.R"

# logs folder relative to where you submit from OR absolute; this is safest:
mkdir -p ${ROOT_DIR}/scripts/PopulationGenetics/logs

echo "Starting rabies ORV choropleth"
date
echo "Running: ${SCRIPT}"

Rscript "${SCRIPT}"

echo "Finished"
date