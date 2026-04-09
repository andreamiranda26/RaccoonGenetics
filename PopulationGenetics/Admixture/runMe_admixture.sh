#!/bin/bash

source $HOME/miniforge3/bin/activate Admixture

ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"

# Define variables
IN_FILE="${ROOT_DIR}/analysis/admixture/admixture.bed"

for K in {1..5}; do
sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=Admx_K${K}
#SBATCH --output=logs/Admixture_K${K}-%j.out
#SBATCH --error=logs/Admixture_K${K}-%j.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --partition=jrw0107_std,general,nova_20
#SBATCH --time=7-00

admixture -j20 -B1000 --cv $IN_FILE $K

mv ${ROOT_DIR}/scripts/PopulationGenetics/Admixture/admixture.${K}* "${ROOT_DIR}/analysis/admixture"
EOF
done