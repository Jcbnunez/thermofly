#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=multiqc 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=2 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu

#--------------------------------------------------------------------------------

# This script will initiate a pipeline which will do some quality QC on the reads.

# Load modules 
# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name multiqc #If you haven't already done so, create and name the environment
source activate multiqc #activate the environment
#conda install -c bioconda multiqc # f you haven't already done so, install the program
conda activate multiqc 

#--------------------------------------------------------------------------------

#Define important file locations

#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

FAST_QC_FOLDER=$WORKING_FOLDER/fastQC

#Name of pipeline
PIPELINE=multiqc

#--------------------------------------------------------------------------------

# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Move to working directory
cd $WORKING_FOLDER

# Generating new folders 
if [ -d "multiqc" ]
then echo "Working multiqc folder exist"; echo "Let's move on"; date
else echo "Working multiqc folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/multiqc; date
fi

#--------------------------------------------------------------------------------

# Run multiqc on the fastqc files

# Move to working directory
cd $WORKING_FOLDER

# Lets do some QC on the reads

# Run multiqc on all of the reads
multiqc $FAST_QC_FOLDER \
-n multiqc_report_all.html \
-o multiqc

# Run multiqc on each of the lanes individually

# Run multiqc on L002 
multiqc $FAST_QC_FOLDER \
-n multiqc_report_L002.html \
--ignore "*L007*" \
--ignore "*L008*" \
-o multiqc

# Run multiqc on L007
multiqc $FAST_QC_FOLDER \
-n multiqc_report_L007.html \
--ignore "*L002*" \
--ignore "*L008*" \
-o multiqc

# Run multiqc on L008 
multiqc $FAST_QC_FOLDER \
-n multiqc_report_L008.html \
--ignore "*L002*" \
--ignore "*L007*" \
-o multiqc

#--------------------------------------------------------------------------------

# This part of the pipeline will produce a notification stating the completion of the script. 

echo "pipeline" ${PIPELINE} $(date)

# Deactivate conda
conda deactivate