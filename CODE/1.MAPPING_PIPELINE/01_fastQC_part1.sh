#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fastQC_part1 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=5 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/fastQC.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu

#--------------------------------------------------------------------------------

# This script will initiate a pipeline which will do some quality QC on the reads.
# Due to the number of fastq files and array limits, the QC is broken up into two arrays 

# Load modules 
# Call fastqc package 
spack load fastqc@0.11.7

#--------------------------------------------------------------------------------

#Define important file locations

#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#Name of pipeline
PIPELINE=fastQC

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.
GUIDE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/Guide_Files/Guide_File_qc_part1.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##               File               Snail_ID  Sample#  Lane# 
## FB1-1_S84_L002_R1_001.fastq.gz    FB1-1     S84     L002
## FB1-1_S84_L002_R2_001.fastq.gz    FB1-1     S84     L002
## FB1-1_S84_L007_R1_001.fastq.gz    FB1-1     S84     L007
## FB1-1_S84_L007_R2_001.fastq.gz    FB1-1     S84     L007
## ...
## HC4-9_S129_L008_R1_001.fastq.gz   HC4-9    S129     L008
## HC4-9_S129_L008_R2_001.fastq.gz   HC4-9    S129     L008

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo ${i}

#--------------------------------------------------------------------------------

# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Move to working directory
cd $WORKING_FOLDER

# Generating new folders 
if [ -d "fastQC" ]
then echo "Working fastQC folder exist"; echo "Let's move on"; date
else echo "Working fastQC folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/fastQC; date
fi

#--------------------------------------------------------------------------------

# Do QC on raw reads with fastqc

# Move to working directory
cd $WORKING_FOLDER

echo -e $i "is now processing"
date

# Lets do some QC on the reads
fastqc $RAW_READS/${i} \
--outdir $WORKING_FOLDER/fastQC

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will produce a notification stating the completion of the script. 

echo ${i} " completed"

echo "pipeline" ${PIPELINE} $(date)