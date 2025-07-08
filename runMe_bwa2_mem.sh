#!/bin/bash

source $HOME/miniforge3/bin/activate Bwa # bwa-mem2 2.2.1

#Hard paths are better to see where everything is to find 
ROOT_DIR="/home/azm0272/Raccoons"

# Define variables
#Grab the genome from Genebank and save it using wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/646/535/GCA_028646535.1_pl-1k/GCA_028646535.1_pl-1k_genomic.fna.gz
#remember to decompress using gunzip so run this to keep the original and write the output separately gunzip -c GCA_028646535.1_pl-1k_genomic.fna.gz > GCA_028646535.1_pl-1k_genomic.fna


SAMPLES=($(ls ${ROOT_DIR}/Willoughby2_Project_003/*fastq.gz | cut -d'/' -f6 | cut -d'_' -f1))

GENOME="GCA_028646535.1_pl-1k_genomic.fna"
GENOME_FILE="${ROOT_DIR}/Raccoon_analysis/data/reference/${GENOME}"
FASTP_DIR="${ROOT_DIR}/Raccoon_analysis/data/fastp"

mkdir -p ${ROOT_DIR}/Raccoon_analysis/data/bwamem2

# Run BWA-MEM2
for SAMPLE in ${SAMPLES[@]}; do

    FASTQ="${FASTP_DIR}/${SAMPLE}.fq.gz"
    BAM="${ROOT_DIR}/Raccoon_analysis/data/bwamem2/${SAMPLE}.bam"

sbatch <<- EOF
#!/bin/bash

#SBATCH --job-name=BWA2-mem_${SAMPLE}
#SBATCH --output=logs/BWA2-mem/BWA2-mem_${SAMPLE}-%j.out          # Output file. %j is replaced with job ID
#SBATCH --cpus-per-task=10
#SBATCH --time=5-00
#SBATCH --mem=40G
#SBATCH -p jrw0107_std,general,nova

bwa-mem2 mem \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina" \
    -t 10 \
    ${GENOME_FILE} \
    ${FASTQ} \
    | samtools sort --threads 10 -o ${BAM} -
EOF
done
