#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Clean_bams

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
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Clean_bams.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will clean the bam files, specifically it will filter, sort, remove duplicates and index. 
# I will also conduct an intermediary QC step with Qualimap. 

# Load modules  
spack load samtools@1.10

PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

#Define important file locations

#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Name of pipeline
PIPELINE=Clean_bams

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
echo ${i}

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

if [ -d "bams_clean" ]
then echo "Working bams_clean folder exist"; echo "Let's move on."; date
else echo "Working bams_clean folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_clean; date
fi

if [ -d "bams_qualimap" ]
then echo "Working bams_qualimap folder exist"; echo "Let's move on."; date
else echo "Working bams_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_qualimap; date
fi

#--------------------------------------------------------------------------------

# Clean the bam files

# These steps will sort the bams, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. 
# Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER

# Filter merged bam files with samtools view and add flags
samtools view \
-b \
-q $QUAL \
-f 0x0002 -F 0x0004 -F 0x0008 \
--threads $CPU  \
$WORKING_FOLDER/bams/${i}.bam \
> $WORKING_FOLDER/bams_clean/${i}.bam
# -q = Skip alignments with MAPQ smaller than $QUAL (40)
# 0x0002 = read mapped in proper pair (0x2)*
# 0x0004 = read unmapped (0x4)
# 0x0008 = mate unmapped (0x8)*

# Sort with picard
# Notice that once a file has been sorted it is added the "srt" suffix
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$WORKING_FOLDER/bams_clean/${i}.bam \
O=$WORKING_FOLDER/bams_clean/${i}.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

# Index with samtools
samtools index $WORKING_FOLDER/bams_clean/${i}.srt.bam

# Remove duplicates with picard
# Notice that once a file has duplicates removed it is added the "rmdp" suffix
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
I=$WORKING_FOLDER/bams_clean/${i}.srt.bam \
O=$WORKING_FOLDER/bams_clean/${i}.srt.rmdp.bam \
M=$WORKING_FOLDER/bams_clean/${i}.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# Lets do QC on the bam file
$qualimap bamqc \
-bam $WORKING_FOLDER/bams_clean/${i}.srt.rmdp.bam \
-outdir $WORKING_FOLDER/bams_qualimap/Qualimap_${i} \
--java-mem-size=$JAVAMEM

# Clean intermediate files
rm $WORKING_FOLDER/bams_clean/${i}.bam
rm $WORKING_FOLDER/bams_clean/${i}.srt.bam

# Housekeeping
mv $WORKING_FOLDER/bams/${i}.dupstat.txt \
$WORKING_FOLDER/mapping_stats

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/Logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)