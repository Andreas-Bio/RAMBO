# RAMBO
Resolving Amplicons in Mixed Samples for Accurate DNA Barcoding with Oxford Nanopore

## Setup and Initialization

### Running the Pipeline

To begin, open a **new R session** (RStudio or terminal) to ensure a clean environment. Then run the main driver script using:

```
source("path/to/denoise.R")
```
or
```
Rscript path/to/denoise.R
```

It is recommended to run the script from a directory where all input and output folders are properly structured, and to avoid workspace contamination from previous runs.
The folder specified in the config file should be empty, except for the designated input files in fastq or fasta format.

### System Requirements

Before starting, make sure the following system-level dependencies are available:

- **GNU Make**: Required for compiling packages that include native code. This is typically available on Unix systems; Windows users must install [RTools](https://cran.r-project.org/bin/windows/Rtools/).
- **MAFFT**: Required for sequence alignment. Ensure that the `mafft` executable is available and set via the `my_mafft` path variable in the script.
- **Internet access** is required to install any missing R packages on first run.

### Directory Setup

Ensure that the following subdirectories exist relative to the main script:
- `scripts/` – contains sourced helper functions and modules

### Script Directory Detection

The pipeline uses the `this.path` package to determine the location of the executing script. This allows file paths to be set relative to the script’s directory:

```
if (!"this.path" %in% installed.packages()) install.packages("this.path")
library(this.path)
```

All necessary R script components are loaded from the `scripts/` directory:

```
source(file = paste0(this.path::here(), "/scripts/startup.R"))
source(file = paste0(this.path::here(), "/scripts/align_fastq.R"))
source(file = paste0(this.path::here(), "/scripts/filter_nanopore_batch.R"))
source(file = paste0(this.path::here(), "/scripts/png_gallery_to_pdf_paginated.R"))
source(file = paste0(this.path::here(), "/scripts/homopolymer.R"))
source(file = paste0(this.path::here(), "/scripts/clustering.R"))
source(file = paste0(this.path::here(), "/scripts/build_tools.R"))
source(file = paste0(this.path::here(), "/scripts/post_merge.R"))
source(file = paste0(this.path::here(), "/scripts/consensus_first_pass.R"))
```

To ensure required build tools are available (especially on Windows), the function `check_build_environment()` verifies that system compilers (e.g., RTools or `make`) are accessible.
If this function fails to detect the RTools installation please provide the correct path manually (see below).

```
check_build_environment()  # Optionally: check_build_environment(rtools_dir = "C:/Rtools43")
```

Required R packages are installed and loaded via a helper function `ensure_required_packages()`, which handles both CRAN and Bioconductor dependencies.
This step can be time-consuming, and if you don't have installation rights to the default R package directory, it may result in an error.
Biostrings-related packages are automatically installed via the Bioconductor package manager.

```
required <- c(
  "ShortRead", "Biostrings", "Rcpp", "ggplot2", "parallel", "doParallel", "foreach",
  "grid", "png", "htmlwidgets", "plotly", "dbscan", "FNN", "igraph", "pwalign",
  "magrittr", "uwot", "Matrix", "dplyr", "stringr"
)

ensure_required_packages(required)
```

This ensures the R session is fully configured before executing any downstream analysis.
