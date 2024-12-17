#!/usr/bin/env bash  
#  
#SBATCH -J HaploCall  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 20:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p bluemoon  
#SBATCH --array=1-15,17-18,21

#--- Notes ---------------------------------------------------------------------
# samples 16, 19, 20 and 22 have too low coverage
#--------------------------------------------------------------------------------
#### Alternative filtering analyses

#### Load Modules
module load singularity
gatk=/netfiles/nunezlab/Shared_Resources/Software/gatk_latest.sif

picard=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip

#### Locations
BAMS_FOLDER=/netfiles/thermofly/bams_clean/
WORKING_FOLDER=/users/j/c/jcnunez/scratch/thermofly/basisetae/mapping
REFERENCE=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked.fa

#### Metadata
meta=/netfiles/thermofly/METADATA/Thermofly_metadata.tsv
SUFFIX=srt.rmdp

#### Java info
JAVAMEM=18G
CPU=$SLURM_CPUS_ON_NODE
echo ${CPU}

#Read Information
Group_library="Thermofly_Don"
Library_Platform="G4"
Group_platform="Thermofly"
#HaploCaller
HET=0.005

#####
i=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo ${i}

###########################################################################
###########################################################################
# Forcing a uniform read group to the joint bam file
###########################################################################
###########################################################################

mkdir $WORKING_FOLDER/RGSM_final_bams

java -jar $picard AddOrReplaceReadGroups \
  I=$BAMS_FOLDER/${i}.$SUFFIX.bam \
  O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
  RGLB=$Group_library \
  RGPL=$Library_Platform \
  RGPU=$Group_platform \
  RGSM=${i}

###########################################################################
###########################################################################
# Index Bam files
###########################################################################
###########################################################################

java -jar $picard BuildBamIndex \
      I=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
      O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bai

###########################################################################
###########################################################################
# Haplotype Calling
###########################################################################
###########################################################################

mkdir $WORKING_FOLDER/haplotype_calling

# NEED TO MAKE DICTIONARY
#java -jar $picard CreateSequenceDictionary \
#      R=$REFERENCE \
#      O=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked.dict

SINGULARITYENV_i=${i} \
SINGULARITYENV_REFERENCE=${REFERENCE} \
SINGULARITYENV_WORKING_FOLDER=${WORKING_FOLDER} \
SINGULARITYENV_JAVAMEM=${JAVAMEM} \
SINGULARITYENV_HET=${HET} \
singularity exec -H ${WORKING_FOLDER} ${gatk} \
gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" \
HaplotypeCaller \
  -R $REFERENCE \
  -I $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
  -O $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf \
  --heterozygosity $HET \
  -ploidy 2 \
  -ERC GVCF

###########################################################################
###########################################################################
# Compress and index with Tabix
###########################################################################
###########################################################################

mkdir 

$bgzip $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf
$tabix $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf.gz

echo "done"

