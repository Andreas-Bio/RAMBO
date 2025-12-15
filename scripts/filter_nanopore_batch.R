validate_nanopore_filter_inputs <- function(fastq_files, n_cores = 1) {
  # Load packages (if not already loaded)
  if (!requireNamespace("ShortRead", quietly = TRUE)) stop("Package 'ShortRead' is required.")
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("Package 'Biostrings' is required.")
  
  # Validate file input
  if (!is.character(fastq_files)) stop("`fastq_files` must be a character vector.")
  fastq_files <- fastq_files[!is.na(fastq_files) & fastq_files != ""]
  if (length(fastq_files) == 0) stop("No FASTQ file paths provided.")
  
  missing_files <- fastq_files[!file.exists(fastq_files)]
  if (length(missing_files) > 0) {
    stop("The following FASTQ files do not exist:\n", paste(missing_files, collapse = "\n"))
  }
  
  # Attempt to read headers from each file to ensure FASTQ format
  valid_fastq <- vapply(fastq_files, function(f) {
    tryCatch({
      obj <- ShortRead::readFastq(f, n = 1)
      is(obj, "ShortReadQ") && length(obj) > 0
    }, error = function(e) FALSE)
  }, logical(1))
  
  if (!all(valid_fastq)) {
    stop("Some files are not valid FASTQ files readable by ShortRead:\n",
         paste(fastq_files[!valid_fastq], collapse = "\n"))
  }
  
  # Validate number of cores
  available_cores <- parallel::detectCores(logical = TRUE)
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1 || n_cores > available_cores) {
    stop(sprintf("`n_cores` must be between 1 and %d (available). Got: %s", available_cores, n_cores))
  }
  
  # Return cleaned inputs
  list(fastq_files = fastq_files, n_cores = as.integer(n_cores))
}

filter_nanopore_batch <- function(
    fastq_files,
    suffix = "_filtered",
    max_removal_fraction = 0.05,
    mad_threshold = 3,
    plot = FALSE,
    n_cores = parallel::detectCores(logical = FALSE)
) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(i = seq_along(fastq_files), .combine = rbind, .packages = c("ShortRead", "ggplot2", "magrittr")) %dopar% {
      tryCatch({
        fq_file <- fastq_files[i] 
        fq_obj <- readFastq(fq_file)
        
        acc <- quality(fq_obj) %>%
          as("matrix") %>%
          { 100 * (1 - 10^(-. / 10)) } %>%
          rowMeans(na.rm = TRUE)
        
        total <- length(acc)
        if (total < 6) {
          return(data.frame(input = basename(fq_file), output = NA, plot = NA,
                            reads_before = total, reads_removed = 0,
                            acc_before = round(mean(acc), 2), acc_after = round(mean(acc), 2),
                            method = "too_few", threshold = NA))
        }
        
        max_rm <- if (total < 10) 1 else floor(total * max_removal_fraction)
        mad_cutoff <- median(acc) - mad_threshold * mad(acc, constant = 1)
        use_q <- sum(acc < mad_cutoff) > max_rm
        threshold <- if (use_q) sort(acc)[max_rm + 1] else mad_cutoff
        keep <- which(acc >= threshold)
        
        if (length(keep) < 5) {
          keep <- order(acc, decreasing = TRUE)[1:5]
          threshold_method <- "forced_min5"
        } else {
          threshold_method <- if (use_q) "quantile_5%" else "mad_based"
        }
        
        out_fq <- sub("(\\.[^.]+)(\\.gz)?$", paste0(suffix, "\\1"), fq_file)
        
        if (file.exists(out_fq)) file.remove(out_fq)
        
        writeFastq(fq_obj[keep], out_fq, compress = FALSE)
        
        out_svg <- sub("\\.fastq(\\.gz)?$", ".png", out_fq, ignore.case = TRUE)
        
        if (plot) {
          plot_obj <- data.frame(acc = acc) %>%
            ggplot(aes(x = acc)) +
            geom_histogram(bins = 60, fill = "skyblue", color = "black") +
            geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
            labs(title = "Read Accuracy", x = "Mean Accuracy (%)", y = "Reads") +
            theme_minimal()
          
          ggsave(
            filename = out_svg,
            plot = plot_obj,
            device = "png",
            width = 6,
            height = 4,
            units = "in"
          )
          
        }
        
        data.frame(
          input = basename(fq_file),
          output = basename(out_fq),
          plot = if (plot) basename(out_svg) else NA,
          reads_before = total,
          reads_removed = total - length(keep),
          acc_before = round(mean(acc), 2),
          acc_after = round(mean(acc[keep]), 2),
          method = threshold_method,
          threshold = round(threshold, 2)
        )
      }, error = function(e) {
        msg <- conditionMessage(e)
        cat(sprintf("ERROR in file [%s]: %s\n", fastq_files[i], msg))
        data.frame(input = basename(fastq_files[i]), output = NA, plot = NA,
                   reads_before = NA, reads_removed = NA,
                   acc_before = NA, acc_after = NA,
                   method = paste("error:", msg), threshold = NA)
      })
    }
  }, finally = {
    stopCluster(cl)  # Always stop the cluster
  })
  
  return(results)
}