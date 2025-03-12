#!/bin/bash -l
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=101-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=Alava.stdout
#SBATCH -o Alava.stdout
#SBTACH -e Alava.stdout
#SBATCH --account jgoudet_pop_fst

module load gcc/11.4.0
module load r-light

#Example
OPTIMA_TYPE="Neutral" 
GENERATIONS=5000
POPULATION_STRUCTURE="SS"
SELECTIVE_OR_NEUTRAL="Neutral_mig0064_5kGen_allreps"
NP=20

Rscript Chapter1/MethodsLocalAdaptation_LogAV/Scripts/Method_Testing_Scripts/MethodsLocalAdaptation_LAVA_QN_test.r \
    $SLURM_ARRAY_TASK_ID \
    $OPTIMA_TYPE \
    $GENERATIONS \
    $POPULATION_STRUCTURE \
    $SELECTIVE_OR_NEUTRAL \
    $NP
