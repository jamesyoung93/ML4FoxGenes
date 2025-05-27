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

1. (Optional) restore R package versions using `renv`:
   ```R
   renv::restore(prompt = FALSE)
   ```
2. Execute the analysis script:
   ```bash
   Rscript analyses/run_analysis.R
   ```

An example SLURM job file `job.slurm` is provided for running on HPC systems.
