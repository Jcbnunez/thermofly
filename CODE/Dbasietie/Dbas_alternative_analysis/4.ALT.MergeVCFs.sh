#!/usr/bin/env bash  
#  
#SBATCH -J Merge  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 12:00:00   
#SBATCH --mem 120G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p bluemoon  

#--- Notes ---------------------------------------------------------------------
# samples 16, 19, 20 and 22 have too low coverage
#--------------------------------------------------------------------------------
#### Alternative filtering analyses

#### Load Modules
picard=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip

#### Locations
WORKING_FOLDER=/gpfs2/scratch/jcnunez/thermofly/basisetae/mapping

#### Namings
PIPELINE=D.basisetae.GATK.pipe

###########################################################################
#Parameters
#Java
JAVAMEM=110G
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
###########################################################################

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Gather VCFs to make a final VCF
###########################################################################
###########################################################################

# order 2L, 2R, 3L, 3R, 4, X, Y

java -Xmx$JAVAMEM \
 -jar $picard GatherVcfs \
  I=$WORKING_FOLDER/contig_1.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_2.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_3.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_4.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_5.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_6.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_7.genotyped.raw.vcf.gz \
  I=$WORKING_FOLDER/contig_8.genotyped.raw.vcf.gz \
  O=$WORKING_FOLDER/$PIPELINE.vcf
	
#bgzip and tabix
	$bgzip $WORKING_FOLDER/$PIPELINE.vcf
	$tabix $WORKING_FOLDER/$PIPELINE.vcf.gz
	
echo "done" $(date)
