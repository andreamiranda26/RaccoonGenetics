#!/bin/bash

#SBATCH --job-name=PlotAdmixture
#SBATCH --output=logs/PlotAdmixture-%j.out          # Output file. %j is replaced with job ID
#SBATCH --error=logs/PlotAdmixture-%j.err           # Error file. %j is replaced with job ID
#SBATCH --cpus-per-task=2
#SBATCH --time=5-00
#SBATCH --mem=10G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate Rmarkdown

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"

cd ${ROOT_DIR}/scripts/PopulationGenetics/Rmd/

Rscript -e "rmarkdown::render('Admixture-plot.Rmd')"
