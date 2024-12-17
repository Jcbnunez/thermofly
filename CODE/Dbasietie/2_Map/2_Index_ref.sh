#!/usr/bin/env bash  
#  
#SBATCH -J Index_ref  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x_%j.out    
#SBATCH -p bluemoon  

#--------------------------------------------------------------------------------

# Load software 
spack load samtools@1.10
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

# Set folders and file locations
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked

#--------------------------------------------------------------------------------

cd /netfiles/thermofly/GENOMES/basisetae

# Index reference

# bwa mem
$bwa index $ref 

# picard
java -jar $PICARD CreateSequenceDictionary \
R=$ref

# samtools
samtools faidx $ref