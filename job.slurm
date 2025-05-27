#!/bin/bash
#SBATCH --job-name=foxgenes
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load R/4.3.2-gcc-11.2.0-mkl
Rscript analyses/run_analysis.R
