#!/bin/bash
#source is a command that says run this script which activates the environment and since you are not activating mamba Fastp 
# the source is a workaround to activate 
source $HOME/miniforge3/bin/activate Fastp

#Hard paths are better to see where everything is to find 
ROOT_DIR="/home/azm0272/Raccoons"

# Define variables
#samples are where the samples are. Here listed every file in the folder that has fastq.gz in the name
#cut by slash and pick file name then cut by underscore and gives me the first part which was sample name 

SAMPLES=($(ls ${ROOT_DIR}/Willoughby2_Project_003/*fastq.gz | cut -d'/' -f6 | cut -d'_' -f1))
#this is where the fastp files are going to be saved 
FASTP_DIR="${ROOT_DIR}/Raccoon_analysis/data/fastp"

mkdir -p ${FASTP_DIR}/{logs,fail}
# Run 
#this is syntax to iterate over a list and the raw is what we name the raw files we got from sequencing 
#EOF is end of file and this limits what is put into sbatch 
for SAMPLE in ${SAMPLES[@]}; do

    RAW="${ROOT_DIR}/Willoughby2_Project_003/${SAMPLE}*.fastq.gz"

sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=Fastp-${SAMPLE}
#SBATCH --output=logs/fastp/Fastp-${SAMPLE}.out          # Output file. %j is replaced with job ID
#SBATCH --cpus-per-task=6
#SBATCH --time=5-00
#SBATCH --mem=5G
#SBATCH -p jrw0107_std,general,nova

fastp --in1 ${RAW} \
    --out1 ${FASTP_DIR}/${SAMPLE}.fq.gz \
    --json ${FASTP_DIR}/logs/${SAMPLE}.fastp.json \
    --html ${FASTP_DIR}/logs/${SAMPLE}.fastp.html \
    --failed_out ${FASTP_DIR}/fail/${SAMPLE}.fail.fq.gz \
    --thread 6 \
    --dont_eval_duplication --cut_right --cut_window_size 10 --cut_mean_quality 20 --trim_poly_g
EOF
done
