#!/bin/bash
#SBATCH --job-name=foxgenes
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load R/4.2.0
Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv'); renv::restore(prompt = FALSE)"
Rscript analyses/run_analysis.R
