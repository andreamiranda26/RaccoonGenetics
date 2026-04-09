#!/usr/bin/env bash
# run_plot_structure.sh
# Simple SBATCH job to render Plot_structure.Rmd and save plots to the analysis/structure/plots folder.

source $HOME/miniforge3/bin/activate Rmarkdown
set -euo pipefail

# ========= CONFIG =========
ROOT_DIR="/scratch/Raccoon_andrea/gabe_trial"
RMD="${ROOT_DIR}/scripts/PopulationGenetics/Structure/Plot_structure.Rmd"
LABELS="${ROOT_DIR}/scripts/PopulationGenetics/labels.csv"
OUTDIR="${ROOT_DIR}/analysis/structure/plots"
LOGDIR="${ROOT_DIR}/analysis/structure/logs"
PARTITIONS="jrw0107_std,nova_20,nova_28,general"
TIME="2-00:00:00"     # run time (D-HH:MM:SS)
# ==========================

mkdir -p "$OUTDIR" "$LOGDIR"

sbatch <<- EOF
#!/usr/bin/env bash
#SBATCH --job-name=Plot_STRUCTURE
#SBATCH --partition=${PARTITIONS}
#SBATCH --time=${TIME}
#SBATCH --ntasks=1
#SBATCH --output=${LOGDIR}/Plot_STRUCTURE-%j.out
#SBATCH --error=${LOGDIR}/Plot_STRUCTURE-%j.err

# If your cluster uses modules:
# module load R
# module load r/4.3.2

# export paths for the R Markdown
export ROOT_DIR="${ROOT_DIR}"
export LABELS_CSV="${LABELS}"

# Render the R Markdown to the plots directory
Rscript -e "dir.create('${OUTDIR}', showWarnings=FALSE, recursive=TRUE);
            rmarkdown::render('${RMD}',
                              output_dir='${OUTDIR}',
                              output_file='Structure_Analysis.html')"
EOF

echo "[submitted] RMarkdown render job submitted!"
echo "Outputs will be in: ${OUTDIR}/Structure_Analysis.html"
