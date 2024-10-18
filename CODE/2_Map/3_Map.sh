#!/usr/bin/env bash  
#  
#SBATCH -J Map_reads  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p bluemoon  
#SBATCH --array=1-22
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Map reads to reference

# Load software  
spack load gcc@9.3.0
spack load samtools@1.10
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2

#--------------------------------------------------------------------------------

echo ${SLURM_ARRAY_TASK_ID}

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked

#--------------------------------------------------------------------------------

# Use metadata file to extract sample names and forward and reverse reads
FIL1=$(cat ${meta} | awk -F '\t' '{print $24}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
FIL2=$(cat ${meta} | awk -F '\t' '{print $25}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

SAMP_NAME=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo $FIL1
echo $FIL2
echo $SAMP_NAME

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6
QUAL=40 

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir sams
mkdir mapping_stats
mkdir bams

#--------------------------------------------------------------------------------

# Map reads to reference
echo "Map" ${SAMP_NAME}

# Map with bwa mem2
$bwa mem -M -t $CPU $ref \
$working_folder/cleaned_reads/${SAMP_NAME}_R1_clean.fq.gz \
$working_folder/cleaned_reads/${SAMP_NAME}_R2_clean.fq.gz \
> $working_folder/sams/${SAMP_NAME}.sam

# Extract summary stats
samtools flagstat --threads $CPU \
$working_folder/sams/${SAMP_NAME}.sam \
> $working_folder/mapping_stats/${SAMP_NAME}.flagstats_raw.sam.txt

#--------------------------------------------------------------------------------

# Build bam files
samtools view -b -q $QUAL --threads $CPU  \
$working_folder/sams/${SAMP_NAME}.sam \
> $working_folder/bams/${SAMP_NAME}.bam

#--------------------------------------------------------------------------------

# Inform completion of pipeline
echo "pipeline completed" $(date)