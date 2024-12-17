#!/usr/bin/env bash  
#  
#SBATCH -J Coverage_bams  
#SBATCH -c 1 
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out  
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate coverage

# Load software 
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
#guide_file=$working_folder/METADATA/Qualimap_bam_list.txt
guide_file=$working_folder/METADATA/Qualimap_bam_list_reduced.txt

#--------------------------------------------------------------------------------

# Parameters for software
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
#mkdir bams_qualimap_multi
mkdir bams_qualimap_multi_reduced

#--------------------------------------------------------------------------------

# Assess quality of bam files
$qualimap multi-bamqc \
-d $guide_file \
-outdir $working_folder/bams_qualimap_multi_reduced \
--java-mem-size=$JAVAMEM
