#!/bin/bash -l
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=SS5kgen4
#SBATCH -o SS5kgen4.stdout
#SBTACH -e SS5kgen4.stderr
#SBATCH --account jgoudet_pop_fst

quantinemo_linux/quantinemo quantiNemo_NEUTRAL_SS_20pop.ini