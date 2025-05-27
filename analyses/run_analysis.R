#!/usr/bin/env Rscript

# Restore packages if renv is used
if (requireNamespace("renv", quietly = TRUE) && file.exists("renv.lock")) {
  renv::restore(prompt = FALSE)
}

rmarkdown::render(
  input = file.path("analyses", "PredictingFOXGenes.Rmd"),
  output_dir = "results"
)
