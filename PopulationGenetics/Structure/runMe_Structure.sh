#!/usr/bin/env bash

source $HOME/miniforge3/bin/activate Structure

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"

# Define variables

MAIN="${ROOT_DIR}/scripts/PopulationGenetics/Structure/mainparams"
EXTRA="${ROOT_DIR}/scripts/PopulationGenetics/Structure/extraparams"

IN_FILE="${ROOT_DIR}/analysis/structure/Raccoon.str"
for k in {1..4}; do
    for r in {1..10}; do
        OUTPUT="${ROOT_DIR}/analysis/structure/Raccoon_K${k}_R${r}.out"
sbatch <<- EOF
#!/usr/bin/env bash
#SBATCH --job-name=Str_K${k}_R${r}             # job name
#SBATCH --ntasks=1                  # number of tasks across all nodes
#SBATCH --partition=jrw0107_std,nova_20,nova_28,general          # name of partition to submit job
#SBATCH --time=10-00:00:00              # Run time (D-HH:MM:SS)
#SBATCH --output=logs/Str_K${k}_R${r}-%j.out          # Output file. %A is replaced with job ID
#SBATCH --error=logs/Str_K${k}_R${r}-%j.err           # Error file. %A is replaced with job ID

structure -K $k -o $OUTPUT -m $MAIN -i $IN_FILE -e $EXTRA
EOF
    done
done