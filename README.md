# ML4FoxGenes

This project explores statistical modeling and machine learning approaches for FOX gene discovery.

## Repository Layout

```
analyses/   - R Markdown and scripts for analysis
data/       - input CSV files
results/    - output predictions and rendered reports
docs/       - additional documentation
```

## Running the Analysis

1. **Restore packages (optional)**

   Use `renv` to install the package versions recorded in `renv.lock`. If you
   do not have write access to the system library, create a personal library and
   point `renv` there:

   ```bash
   mkdir -p ~/R/4.3.2-library
   Rscript -e "install.packages('renv', repos='https://cloud.r-project.org', lib='~/R/4.3.2-library')"
   Rscript -e ".libPaths('~/R/4.3.2-library'); renv::restore(prompt = FALSE)"
   ```

2. **Run the analysis**

   The script `analyses/run_analysis.R` renders the R Markdown file when Pandoc
   is available. On systems without Pandoc it automatically falls back to running
   a plain R script extracted from the notebook.

   ```bash
   Rscript analyses/run_analysis.R
   ```

An example SLURM job file `job.slurm` is provided for running on HPC systems.
