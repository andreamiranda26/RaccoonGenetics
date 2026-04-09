#!/bin/bash

#SBATCH --job-name=IBD_raccoon
#SBATCH --output=logs/IBD-%j.out
#SBATCH --error=logs/IBD-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00
#SBATCH --mem=16G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate Rmarkdown

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"

cd ${ROOT_DIR}/scripts/PopulationGenetics/Rmd/


Rscript -e "rmarkdown::render('IBD_mantel.Rmd')"

