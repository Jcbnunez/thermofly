#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=qualimap_multi_bamqc

# Specify partition
#SBATCH --partition=bigmem

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=500G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate average coverage 

#Load modules 
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#--------------------------------------------------------------------------------
# Define parameters
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Move to working folder
cd $WORKING_FOLDER

# Generate output folder 
if [ -d "qualimap_multi_bamqc" ]
then
echo "Working qualimap_multi_bamqc folder exist"
echo "lets move on"
else 
echo "Folder doesnt exist. Lets fix that"
mkdir qualimap_multi_bamqc
fi

OUTPUT=$WORKING_FOLDER/qualimap_multi_bamqc

#--------------------------------------------------------------------------------

## Read guide files
# This is a file with the name all the samples to be processed and the path to each Qualimap.

GUIDE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/Guide_Files/Qualimap_bam_list.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "bams_merged_qualimap_multi" ]
then echo "Working bams_merged_qualimap_multi folder exist"; echo "Let's move on."; date
else echo "Working bams_merged_qualimap_multi folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_merged_qualimap_multi; date
fi

#--------------------------------------------------------------------------------

# Assess quality of all bam files
$qualimap multi-bamqc \
-d $GUIDE_FILE \
-outdir $WORKING_FOLDER/bams_merged_qualimap_multi \
--java-mem-size=$JAVAMEM