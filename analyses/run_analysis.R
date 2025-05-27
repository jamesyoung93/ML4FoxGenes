#!/usr/bin/env Rscript

# Restore packages if renv is used
if (requireNamespace("renv", quietly = TRUE) && file.exists("renv.lock")) {
  renv::restore(prompt = FALSE)
}

input_rmd <- file.path("analyses", "PredictingFOXGenes.Rmd")
script_r <- file.path("analyses", "PredictingFOXGenes.R")
output_dir <- "results"

run_as_script <- function() {
  if (!file.exists(script_r)) {
    message("Extracting R script from Rmd ...")
    knitr::purl(input_rmd, output = script_r)
  }
  message("Running analysis as plain R script")
  source(script_r)
}

if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
  rmarkdown::render(input = input_rmd, output_dir = output_dir)
} else {
  run_as_script()
}
