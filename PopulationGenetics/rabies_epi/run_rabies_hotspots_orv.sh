#!/bin/bash
#SBATCH --job-name=RabiesORV
#SBATCH --output=logs/RabiesORV-%j.out
#SBATCH --error=logs/RabiesORV-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
cd ${ROOT_DIR}/scripts/PopulationGenetics/rabies_epi/
mkdir -p logs

date
Rscript rabies_hotspots_orv_centroid.R
date