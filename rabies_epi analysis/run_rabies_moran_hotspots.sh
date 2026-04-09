#!/bin/bash
#SBATCH --job-name=RabiesMoran
#SBATCH --output=logs/RabiesMoran-%j.out
#SBATCH --error=logs/RabiesMoran-%j.err
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

echo "Starting rabies Moran's I hotspot analysis"
date

Rscript rabies_moran_hotspots.R

echo "Finished"
date