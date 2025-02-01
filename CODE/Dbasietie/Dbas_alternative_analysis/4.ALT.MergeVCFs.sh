#!/usr/bin/env bash  
#  
#SBATCH -J Merge  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 12:00:00   
#SBATCH --mem 120G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  

#--- Notes ---------------------------------------------------------------------
# samples 16, 19, 20 and 22 have too low coverage
#--------------------------------------------------------------------------------
#### Alternative filtering analyses

#### Load Modules
picard=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
tabix=/netfiles/nunezlab/Shared_Resources/Software/htslib-1.21/tabix
bgzip=/netfiles/nunezlab/Shared_Resources/Software/htslib-1.21/bgzip

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
I=$WORKING_FOLDER/JAWNLB010000001.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000002.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000003.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000004.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000005.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000006.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000007.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000008.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000009.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000010.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000011.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000012.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000013.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000014.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000015.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000016.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000017.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000018.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000019.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000020.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000021.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000022.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000023.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000024.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000025.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000026.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000027.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000028.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000029.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000030.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000031.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000032.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000033.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000034.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000035.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000036.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000037.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000038.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000039.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000040.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000041.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000042.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000043.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000044.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000045.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000046.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000047.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000048.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000049.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000050.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000051.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000052.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000053.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000054.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000055.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000056.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000057.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000058.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000059.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000060.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000061.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000062.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000063.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000064.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000065.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000066.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000067.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000068.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000069.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000070.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000071.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000072.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000073.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000074.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000075.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000076.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000077.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000078.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000079.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000080.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000081.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000082.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000083.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000084.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000085.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000086.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000087.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000088.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000089.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000090.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000091.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000092.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000093.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000094.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000095.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000096.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000097.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000098.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000099.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000100.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000101.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000102.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000103.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000104.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000105.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000106.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000107.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000108.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000109.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000110.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000111.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000112.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000113.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000114.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000115.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000116.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000117.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000118.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000119.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000120.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000121.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000122.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000123.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000124.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000125.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000126.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000127.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000128.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000129.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000130.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000131.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000132.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000133.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000134.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000135.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000136.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000137.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000138.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000139.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000140.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000141.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000142.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000143.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000144.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000145.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000146.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000147.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000148.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000149.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000150.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000151.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000152.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000153.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000154.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000155.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000156.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000157.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000158.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000159.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000160.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000161.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000162.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000163.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000164.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000165.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000166.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000167.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000168.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000169.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000170.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000171.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000172.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000173.1.genotyped.raw.vcf.gz \
I=$WORKING_FOLDER/JAWNLB010000174.1.genotyped.raw.vcf.gz \
O=$WORKING_FOLDER/$PIPELINE.vcf
	
#bgzip and tabix
	$bgzip $WORKING_FOLDER/$PIPELINE.vcf
	$tabix $WORKING_FOLDER/$PIPELINE.vcf.gz
	
echo "done" $(date)
