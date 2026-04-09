#!/bin/bash
#SBATCH --job-name=RabiesEpi
#SBATCH --output=logs/RabiesEpi-%j.out
#SBATCH --error=logs/RabiesEpi-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
cd ${ROOT_DIR}/scripts/PopulationGenetics/rabies_epi/

mkdir -p logs

echo "Starting rabies spatial epi + suitability"
date

Rscript rabies_suitability_maps.R

echo "Finished"
date