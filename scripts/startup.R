ensure_required_packages <- function(required_packages) {
  # Bioconductor packages (must be installed with BiocManager)
  bioc_packages <- c("ShortRead", "Biostrings", "pwalign")

  # Helper: Unload any loaded packages
  unload_packages <- function(pkgs) {
    for (pkg in pkgs) {
      if (paste("package", pkg, sep = ":") %in% search()) {
        try(detach(paste("package", pkg, sep = ":"), unload = TRUE, character.only = TRUE), silent = TRUE)
      }
    }
  }

  # Helper: Install a package with retry
  safe_install <- function(pkg, use_bioc = FALSE) {
    max_tries <- 3
    attempt <- 1
    success <- FALSE
    while (attempt <= max_tries && !success) {
      tryCatch({
        if (use_bioc) {
          BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
        } else {
          install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        }
        success <- TRUE
      }, error = function(e) {
        message(sprintf("Installation failed for %s. Retrying (%d/%d)...", pkg, attempt, max_tries))
        attempt <<- attempt + 1
        Sys.sleep(2)
      })
    }
    if (!success) stop(sprintf("Failed to install package: %s after %d attempts", pkg, max_tries))
  }

  # Unload first
  unload_packages(required_packages)

  # Ensure BiocManager is available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
  }

  # Install missing packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      is_bioc <- pkg %in% bioc_packages
      safe_install(pkg, use_bioc = is_bioc)
    }
  }

  # Load everything
  invisible(lapply(required_packages, function(pkg) library(pkg, character.only = TRUE)))
}