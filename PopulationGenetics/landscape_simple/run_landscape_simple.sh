#!/bin/bash
#SBATCH --job-name=LandCoverMantel
#SBATCH --output=logs/LandCoverMantel-%j.out
#SBATCH --error=logs/LandCoverMantel-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
cd ${ROOT_DIR}/scripts/PopulationGenetics/landscape_simple/

mkdir -p logs

echo "Starting simple landscape partial Mantel"
date

Rscript landscape_partial_mantel_county.R

echo "Finished"
date