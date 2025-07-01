#!/usr/bin/env bash  
#  
#SBATCH -J run_glm  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 9:00:00   
#SBATCH --mem 50G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-885

module load Rtidyverse  

n="${SLURM_ARRAY_TASK_ID}"

echo $n

Rscript \
--vanilla \
GLM_PCA_pop_based.R $n

date
echo "done"

