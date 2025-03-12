#!/bin/bash -l
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=ComparedSSsmallpop
#SBATCH -o ComparedSSsmallpop.stdout
#SBTACH -e ComparedSSsmallpop.stderr
#SBATCH --account jgoudet_pop_fst

module load r-light
Rscript Chapter1/MethodsLocalAdaptation_AllMethodsComparison.r $SLURM_ARRAY_TASK_ID