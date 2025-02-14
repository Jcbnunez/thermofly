#!/usr/bin/env bash  
#  
#SBATCH -J Genotype  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 20:00:00   
#SBATCH --mem 160G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-175

###########################################################################
#Parameters
#Java
JAVAMEM=150G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

#Load Modules
#gatk=/netfiles/nunezlab/Shared_Resources/Software/gatk-4.6.0.0/gatk
module load gatk/4.6.1.0

# User defined inputs -- this represents the name of the samples
intervals=/gpfs2/scratch/jcnunez/thermofly/basisetae/Dbas.Intervals.txt

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/users/j/c/jcnunez/scratch/thermofly/basisetae/mapping
REFERENCE=/netfiles/thermofly/GENOMES/basisetae/GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

# Identify the Genome database to genotyoe
GenomeDB_path=`echo $WORKING_FOLDER/db_Drosophila_basisetae/DB_${i}`

echo "now processing DB" ${i} $(date)

###########################################################################
###########################################################################
# Genotype call the samples in the DBI merged set
###########################################################################
###########################################################################

$GATK \
gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
    GenotypeGVCFs \
   -R $REFERENCE \
   -V gendb://$GenomeDB_path \
   -O $WORKING_FOLDER/${i}.genotyped.raw.vcf.gz

echo ${i} "done" $(date)
