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

## Running on an HPC Cluster

Submit the included job script using `sbatch`:

```bash
sbatch job.slurm
```

The script loads R via the module system and restores the package environment
using `renv` before executing `analyses/run_analysis.R`. Ensure your cluster has
network access for package installation or pre-install the required packages in
your user library.
