#!/bin/bash -l
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2
#SBATCH --mem=30G
#SBATCH --job-name=DF_ss_100traits
#SBATCH -o DF_ss_100traits.stdout
#SBTACH -e DF_ss_100traits.stderr
#SBATCH --account jgoudet_pop_fst

module load gcc/11.4.0
module load r
#Example
OPTIMA_TYPE="Neutral" 
GENERATIONS=5000
POPULATION_STRUCTURE="SS_20pop"
SELECTIVE_OR_NEUTRAL="Neutral_5kGen"
NP=20

Rscript MethodsLocalAdaptation_Driftsel_QN_test.r \
    $SLURM_ARRAY_TASK_ID \
    $OPTIMA_TYPE \
    $GENERATIONS \
    $POPULATION_STRUCTURE \
    $SELECTIVE_OR_NEUTRAL \
    $NP
