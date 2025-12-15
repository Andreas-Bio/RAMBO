# ---------------------------------------------------------------------

map_quality <- function(aln_seqs, orig_quals, gap_char) {
  seq_split <- strsplit(aln_seqs, split = "")
  qual_split <- strsplit(orig_quals, split = "")

  mapply(function(seq_vec, qual_vec) {
    gap_mask <- seq_vec == "-"
    out <- rep.int(gap_char, length(seq_vec))
    out[!gap_mask] <- head(qual_vec, sum(!gap_mask))
    paste(out, collapse = "")
  }, seq_split, qual_split, SIMPLIFY = TRUE, USE.NAMES = FALSE)
}


# ---------------------------------------------------------------------
# Align a single FASTQ file using MAFFT
## Patched align_fastq_with_mafft with consistent IDs --------------------
align_fastq_with_mafft <- function(
    fastq_file, path_mafft = "mafft", output_suffix = "_aln",
    param_retree = 2, param_maxiterate = 0, param_threads = 1,
    param_threads_tb = param_threads, param_op = 1, param_ep = 0.12,
    param_kimura = 1, use_reorder = TRUE, additional_mafft_args = character(),
    dir_temp = tempdir(), verbose = TRUE, quiet = TRUE,
    gap_quality_char = "!", overwrite = FALSE, min_avg_qual = NULL
) {
  if (file.info(fastq_file)$size < 100) {
    if (verbose) message("Skipping empty or too-small file: ", fastq_file)
    return(invisible(NULL))
  }
  
  final_fastq_out <- gsub("\\.fastq$", paste0(output_suffix, ".fastq"),
                          fastq_file, ignore.case = TRUE)
  if (file.exists(final_fastq_out) && !overwrite) {
    if (verbose) message("Skipping ", final_fastq_out, " (already exists)")
    return(invisible(NULL))
  }
  
  
  Rcpp::cppFunction('
NumericVector ultra_avg_qscore(CharacterVector qvec) {
  const double min_err = 1e-10;
  double err_lookup[94];
  for (int i = 0; i < 94; ++i) {
    err_lookup[i] = pow(10.0, -i / 10.0);
  }
  R_xlen_t n = qvec.size();
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (qvec[i] == NA_STRING) {
      out[i] = NA_REAL;
      continue;
    }
    std::string s = Rcpp::as<std::string>(qvec[i]);
    double sum_p = 0.0;
    size_t len = s.size();
    bool valid = true;
    for (size_t j = 0; j < len; ++j) {
      int phred = static_cast<unsigned char>(s[j]) - 33;
      if (phred < 0 || phred >= 94) {
        out[i] = NA_REAL;
        valid = false;
        break;
      }
      sum_p += err_lookup[phred];
    }
    if (valid) {
      double mean_err = sum_p / len;
      double q = -10.0 * log10(mean_err > min_err ? mean_err : min_err);
      out[i] = q;
    }
  }
  return out;
}
')
  
  
  fq <- readFastq(fastq_file)
  
  ## optional quality filter -------------------------------------------------
  if (!is.null(min_avg_qual)) {
    quals0_str <- as.character(quality(quality(fq)))
    avg_q <- ultra_avg_qscore(quals0_str)
    keep <- avg_q >= min_avg_qual
    if (verbose) message(sum(!keep), " reads filtered due to low quality")
    if (sum(keep) == 0) {
      if (verbose) message("All reads filtered due to low quality: ", fastq_file)
      return(invisible(NULL))
    }
    if (sum(keep) <= 10) {
      if (verbose) message("Too few reads after filtering (<=10): ", fastq_file)
      return(invisible(NULL))
    }
    fq <- fq[keep]
  }
  
  ## build trimmed ids, sequences, and qualities -----------------------------
  ids_raw  <- as.character(ShortRead::id(fq))
  ids_trim <- sub(" .*", "", ids_raw)
  if (anyDuplicated(ids_trim)) {
    stop("Non-unique trimmed ids; cannot map qualities reliably")
  }
  
  seqs <- ShortRead::sread(fq)
  names(seqs) <- ids_trim
  
  quals     <- quality(fq)
  quals_str <- as.character(quals@quality)
  names(quals_str) <- ids_trim
  
  ## summary stats before alignment -----------------------------------------
  avg_len_before <- mean(width(seqs))
  min_len_before <- min(width(seqs))
  avg_qual_before <- mean(ultra_avg_qscore(as.character(quality(quality(fq)))),
                          na.rm = TRUE)
  seq_count <- length(seqs)
  
  ## write temporary FASTA for MAFFT ----------------------------------------
  temp_fasta_in  <- tempfile(fileext = ".fasta")
  temp_fasta_out <- tempfile(fileext = ".fasta")
  writeXStringSet(seqs, temp_fasta_in)
  
  mafft_args <- c(
    "--retree", param_retree,
    "--maxiterate", param_maxiterate,
    "--thread", param_threads,
    "--threadtb", param_threads_tb,
    "--op", param_op,
    "--ep", param_ep,
    "--kimura", param_kimura,
    if (use_reorder) "--reorder",
    additional_mafft_args,
    temp_fasta_in
  )
  
  system2(path_mafft, args = mafft_args,
          stdout = temp_fasta_out,
          stderr = if (quiet) NULL else "")
  
  if (!file.exists(temp_fasta_out)) {
    if (verbose) message("MAFFT failed to produce output: ", temp_fasta_out)
    unlink(temp_fasta_in)
    return(invisible(NULL))
  }
  
  aln <- readDNAStringSet(temp_fasta_out)
  aln_ids_raw  <- names(aln)
  aln_ids_trim <- sub(" .*", "", aln_ids_raw)
  names(aln)   <- aln_ids_trim
  
  if (length(setdiff(aln_ids_trim, ids_trim)) > 0L) {
    warning("Some aligned ids do not match trimmed FASTQ ids; dropping them")
  }
  
  ## subset qualities to aligned ids ----------------------------------------
  qual_sub <- quals_str[aln_ids_trim]
  if (any(is.na(qual_sub))) {
    bad <- aln_ids_trim[is.na(qual_sub)]
    stop(sprintf("Missing quality strings for %d aligned sequences, e.g. %s",
                 length(bad), bad[1L]))
  }
  
  ## map qualities into gapped alignment, using existing map_quality() -------
  aln_vec   <- as.character(aln)
  aln_quals <- map_quality(aln_vec, qual_sub, gap_char = gap_quality_char)
  
  ## final length sanity check ----------------------------------------------
  if (!all(width(aln) == nchar(aln_quals))) {
    stop("sread and quality widths differ after mapping; aborting")
  }
  
  ## recover full original ids for output -----------------------------------
  id_out_full <- ids_raw[match(aln_ids_trim, ids_trim)]
  
  aligned_fq <- ShortReadQ(
    sread   = aln,
    quality = PhredQuality(BStringSet(aln_quals)),
    id      = BStringSet(id_out_full)
  )
  
  if (file.exists(final_fastq_out) && overwrite) file.remove(final_fastq_out)
  writeFastq(aligned_fq, final_fastq_out, compress = FALSE)
  
  if (verbose) message("Aligned: ", basename(fastq_file))
  
  gap_perc <- sum(letterFrequency(aln, "-")) / sum(width(aln))
  
  unlink(temp_fasta_in)
  unlink(temp_fasta_out)
  
  cbind.data.frame(
    seq_count                 = seq_count,
    len_after_aln             = width(aln)[1],
    avg_len_before_aln        = round(avg_len_before),
    min_len_before_aln        = min_len_before,
    gap_perc_after_aln        = round(gap_perc * 100),
    avg_qual_score_before_aln = round(avg_qual_before)
  )
}
