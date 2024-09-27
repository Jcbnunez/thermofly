#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Merge_bams

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Submit job array
#SBATCH --array=1-192%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Merge_bams.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will gather all sample data across the many lanes of sequencing

# Load modules  
spack load samtools@1.10

PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Path to the directory with the bams (filtered, sorted and duplicates removed). 
BAMS=$WORKING_FOLDER/bams_clean

#Name of pipeline
PIPELINE=Merge_bams

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

GUIDE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/Guide_Files/Guide_File_merge.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Snail_ID  Sample#   Merged_name 1    Merged_name 2    Merged_name 3   Merged_bam_name
##  FB1-1     S84     FB1-1_S84_L002   FB1-1_S84_L007   FB1-1_S84_L008     FB1-1_S84
##  FB1-2     S173    FB1-2_S173_L002  FB1-2_S173_L007  FB1-2_S173_L008    FB1-2_S173
##  FB1-5     S109    FB1-5_S109_L002  FB1-5_S109_L007  FB1-5_S109_L008    FB1-5_S109
##  ...
##  MP9-9     S191    MP9-9_S191_L002  MP9-9_S191_L007  MP9-9_S191_L008    MP9-9_S191
##  MP9-10    S26     MP9-10_S26_L002  MP9-10_S26_L007  MP9-10_S26_L008    MP9-10_S26


#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i

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

if [ -d "bams_merged" ]
then echo "Working bams_merged folder exist"; echo "Let's move on."; date
else echo "Working bams_merged folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_merged; date
fi

if [ -d "bams_merged_qualimap" ]
then echo "Working bams_merged_qualimap folder exist"; echo "Let's move on."; date
else echo "Working bams_merged_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_merged_qualimap; date
fi

#--------------------------------------------------------------------------------

# Start pipeline

# Here I will merge the bam outputs for the multiple lanes of sequencing. These will be named 'Lanes merged'

echo "I will merge these files" $BAMS/${i}_*.srt.rmdp.bam

#Make temporary linefile with list of input BAM files
ls $BAMS/${i}_*.srt.rmdp.bam > ${i}.guide.txt

samtools merge \
-b ${i}.guide.txt \
$WORKING_FOLDER/bams_merged/${i}.Lanes_merged.bam

#remove the temporary guide file
rm ${i}.guide.txt

# Assess quality of final file
$qualimap bamqc \
-bam $WORKING_FOLDER/bams_merged/${i}.Lanes_merged.bam \
-outdir $WORKING_FOLDER/bams_merged_qualimap/Qualimap_LaneMerged_${i} \
--java-mem-size=$JAVAMEM

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/Logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)

