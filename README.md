# RAMBO
Resolving Amplicons in Mixed samples for Barcoding with ONT

Setup and Initialization
This section prepares the environment for running the pipeline. It ensures all required packages are installed, loads core dependencies, sources supporting scripts, and validates the system's build environment.

Script Directory Detection
The script uses the this.path package to dynamically determine the file location of the executing script. This allows all paths to be constructed relative to the script’s directory, avoiding the need for hardcoded file paths.

r
Copy
Edit
if (!"this.path" %in% installed.packages()) install.packages("this.path")
library(this.path)
Sourcing Modular Functions
All functional modules required for the pipeline are sourced from the scripts/ directory. These include filtering, alignment, clustering, consensus generation, and utility functions.

r
Copy
Edit
source(file = paste0(this.path::here(), "/scripts/startup.R"))
source(file = paste0(this.path::here(), "/scripts/align_fastq.R"))
source(file = paste0(this.path::here(), "/scripts/filter_nanopore_batch.R"))
source(file = paste0(this.path::here(), "/scripts/png_gallery_to_pdf_paginated.R"))
source(file = paste0(this.path::here(), "/scripts/homopolymer.R"))
source(file = paste0(this.path::here(), "/scripts/clustering.R"))
source(file = paste0(this.path::here(), "/scripts/build_tools.R"))
source(file = paste0(this.path::here(), "/scripts/post_merge.R"))
source(file = paste0(this.path::here(), "/scripts/consensus_first_pass.R"))
Build Environment Check
The function check_build_environment() verifies that required build tools (e.g. make, compilers) are available. On Windows, it ensures that RTools is installed and properly configured. On Unix-based systems, it checks for the presence of a make utility. The function can optionally take a custom RTools path as an argument.

r
Copy
Edit
check_build_environment()  # Optionally: check_build_environment(rtools_dir = "C:/Rtools43")
Package Installation and Loading
All required CRAN and Bioconductor packages are checked, and if missing, are installed using a helper function ensure_required_packages(). This function should be defined in startup.R or a separate helper script.

r
Copy
Edit
required <- c(
  "ShortRead", "Biostrings", "Rcpp", "ggplot2", "parallel", "doParallel", "foreach",
  "grid", "png", "htmlwidgets", "plotly", "dbscan", "FNN", "igraph", "pwalign",
  "magrittr", "uwot", "Matrix", "dplyr", "stringr"
)

ensure_required_packages(required)
This ensures that the R session has all dependencies installed before the pipeline begins processing data.

