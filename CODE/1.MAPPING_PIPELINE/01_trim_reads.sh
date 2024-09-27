#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Trim_reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=6

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Trim_reads.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will trim the reads.

# Load modules  
fastp=/gpfs1/home/e/l/elongman/software/fastp

#--------------------------------------------------------------------------------

#Define important file locations

#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#Name of pipeline
PIPELINE=Trim_reads

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

GUIDE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/Guide_Files/Guide_File_trim_map.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##               File1                             File2              Snail_ID  Sample#  Lane#    Paired_name    Bam_lanes_merged_name
## FB1-1_S84_L002_R1_001.fastq.gz    FB1-1_S84_L002_R2_001.fastq.gz    FB1-1     S84     L002    FB1-1_S84_L002       FB1-1_S84
## FB1-1_S84_L007_R1_001.fastq.gz    FB1-1_S84_L007_R2_001.fastq.gz    FB1-1     S84     L007    FB1-1_S84_L007       FB1-1_S84
## FB1-1_S84_L008_R1_001.fastq.gz    FB1-1_S84_L008_R2_001.fastq.gz    FB1-1     S84     L008    FB1-1_S84_L008       FB1-1_S84
## FB1-2_S173_L002_R1_001.fastq.gz   FB1-2_S173_L002_R2_001.fastq.gz   FB1-2     S173    L002    FB1-2_S173_L002      FB1-2_S173
## ...
## MP9-10_S26_L007_R1_001.fastq.gz   MP9-10_S26_L007_R2_001.fastq.gz   MP9-10    S26     L007    MP9-10_S26_L007      MP9-10_S26
## MP9-10_S26_L008_R1_001.fastq.gz   MP9-10_S26_L008_R2_001.fastq.gz   MP9-10    S26     L008    MP9-10_S26_L008      MP9-10_S26

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
read1=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
read2=`awk -F "\t" '{print $2}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo ${i} "+" ${read1} "+" ${read2}

#--------------------------------------------------------------------------------

# This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER
# Move to logs direcotry
cd Logs

echo $PIPELINE

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/Logs/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "trimmed_reads" ]
then echo "Working trimmed_reads folder exist"; echo "Let's move on."; date
else echo "Working trimmed_reads folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/trimmed_reads; date
fi

if [ -d "fastp_reports" ]
then echo "Working fastp_reports folder exist"; echo "Let's move on."; date
else echo "Working fastp_reports folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastp_reports; date
fi

#--------------------------------------------------------------------------------

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER

# This script uses an array and matches left and right reads and cleans the raw data using fastp according to parameters set below
# Decide on the trimming parameters based on fastQC step done before this script.

echo ${i} "Trimming reads"

# Call fastp and do some light trimming
$fastp \
-i $RAW_READS/${read1} \
-I $RAW_READS/${read2} \
-o $WORKING_FOLDER/trimmed_reads/${i}_R1_clean.fq.gz \
-O $WORKING_FOLDER/trimmed_reads/${i}_R2_clean.fq.gz \
--detect_adapter_for_pe \
--trim_front1 12 \
--trim_poly_g \
--thread 6 \
--cut_right \
--cut_right_window_size 6 \
--qualified_quality_phred 20 \
--html $WORKING_FOLDER/fastp_reports/${i}_clean.html \
--json $WORKING_FOLDER/fastp_reports/${i}_clean.json

# i = read 1
# I = read 2
# o & O = outputs
# outputs .json and html for qc

# Adapter trimming is enabled by default using overlap analysis, but for PE data you can also specify adapter sequence auto-detection by specifying --detect_adapter_for_pe
# For read 1 of PE data, the front trimming settings are --trim_front1
# Detect and trim polyG in read tails using --trim_poly_g
# Per read cutting by quality options:
### --cut_right: Move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, then stop
### --cut_right_window_size: The window size option for cut_right
### --qualified_quality_phred: The quality vlue that a base is qualified 

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/Logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)