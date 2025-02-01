#!/usr/bin/env bash  
#  
#SBATCH -J HaploCall  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 20:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-15,17-18,21

#--- Notes ---------------------------------------------------------------------
# samples 16, 19, 20 and 22 have too low coverage
#--------------------------------------------------------------------------------
#### Alternative filtering analyses

#### Load Modules
#gatk=/netfiles/nunezlab/Shared_Resources/Software/gatk_latest.sif
module load gatk/4.6.1.0
module load gcc/13.3.0-xp3epyt picard/3.1.1-otrgwkh
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib-1.21/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib-1.21/bgzip

#### Locations
BAMS_FOLDER=/netfiles/thermofly/bams_clean_DPrice_Dbas_Jan2025
WORKING_FOLDER=/users/j/c/jcnunez/scratch/thermofly/basisetae/mapping
REFERENCE=/netfiles/thermofly/GENOMES/basisetae/GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa

#### Metadata
meta=/netfiles/thermofly/METADATA/Thermofly_metadata.vNov11.2024.tsv
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

picard AddOrReplaceReadGroups \
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

picard BuildBamIndex \
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

$GATK \
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

$bgzip $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf
$tabix $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf.gz

echo "done"

