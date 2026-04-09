#!/bin/bash
#SBATCH --job-name=mapRabies
#SBATCH --output=logs/mapRabies-%j.out
#SBATCH --error=logs/mapRabies-%j.err
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:20
#SBATCH --mem=4G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate Rmarkdown

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
cd ${ROOT_DIR}/scripts/PopulationGenetics || exit 1

Rscript map_rabies_orv_baseMaps.R