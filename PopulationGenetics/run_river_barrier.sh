#!/bin/bash
#SBATCH --job-name=riverBarrier
#SBATCH --output=logs/riverBarrier-%j.out
#SBATCH --error=logs/riverBarrier-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH -p jrw0107_std,general,nova_20,nova_28
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=azm0272@auburn.edu

source $HOME/miniforge3/bin/activate raccoon_spatial

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
SCRIPT="${ROOT_DIR}/scripts/PopulationGenetics/river_barrier_test.R"

mkdir -p ${ROOT_DIR}/scripts/PopulationGenetics/logs

echo "Starting river barrier test"
date

Rscript ${SCRIPT}

echo "Finished"
date