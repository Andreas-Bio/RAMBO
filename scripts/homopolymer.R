filter_homopolymer_alignment <- function(
    fastq_path,
    fastq_dir             = "./",
    min_len               = 6,
    gap_fraction_cutoff   = 0.2,
    homopolymer_threshold = 0.5,
    di_rep                = 4L,
    tri_rep               = 4L,
    use_di_tri            = FALSE,
    output_suffix         = "_nohomo"
){
  # C++ Homopolymer matrix (mono runs)
  if (!exists("get_hpmatrix_cpp", mode = "function")) {
    Rcpp::cppFunction('
      IntegerMatrix get_hpmatrix_cpp(CharacterMatrix seq_mat, int min_len = 6) {
        int nrow = seq_mat.nrow();
        int ncol = seq_mat.ncol();
        IntegerMatrix hpmat(nrow, ncol);
        for (int i = 0; i < nrow; ++i) {
          CharacterVector row = seq_mat(i, _);
          std::string run_base = "";
          std::vector<int> run_indices;
          int base_count = 0;
          for (int j = 0; j < ncol; ++j) {
            std::string base = as<std::string>(row[j]);
            if (base == "A" || base == "C" || base == "G" || base == "T") {
              if (run_base == "") {
                run_base = base;
                run_indices = {j};
                base_count = 1;
              } else if (base == run_base) {
                run_indices.push_back(j);
                base_count++;
              } else {
                if (base_count >= min_len) {
                  for (int idx : run_indices) {
                    hpmat(i, idx) = 1;
                  }
                }
                run_base = base;
                run_indices = {j};
                base_count = 1;
              }
            } else if (base == "-" && run_base != "") {
              run_indices.push_back(j);  // gap within run
            } else {
              if (base_count >= min_len) {
                for (int idx : run_indices) {
                  hpmat(i, idx) = 1;
                }
              }
              run_base.clear();
              run_indices.clear();
              base_count = 0;
            }
          }
          if (base_count >= min_len) {
            for (int idx : run_indices) {
              hpmat(i, idx) = 1;
            }
          }
        }
        return hpmat;
      }
    ')
  }
  
  fq    <- ShortRead::readFastq(fastq_path)
  seqs  <- toupper(as.character(ShortRead::sread(fq)))
  stopifnot(length(unique(nchar(seqs))) == 1)
  
  seq_mat   <- do.call(rbind, strsplit(seqs, split = ""))
  phred_mat <- as(Biostrings::quality(fq), "matrix")
  
  # 1) homopolymer columns (mono)
  hpmat      <- get_hpmatrix_cpp(seq_mat, min_len = min_len)
  hpcol_mask <- colSums(hpmat) / nrow(hpmat) >= homopolymer_threshold
  
  # 2) di and tri repeats (optional)
  rep_mask <- rep(FALSE, ncol(seq_mat))
  
  if (use_di_tri) {
    bases <- c("A","C","G","T")
    di    <- as.vector(outer(bases, bases, paste0))                       # AA, AC, ..., TT
    tri   <- as.vector(outer(outer(bases, bases, paste0), bases, paste0)) # AAA, AAC, ..., TTT
    
    dinu   <- vapply(di,  function(m) paste0(rep(m,  di_rep),  collapse = ""), "")
    triu   <- vapply(tri, function(m) paste0(rep(m, tri_rep), collapse = ""),  "")
    motifs <- c(dinu, triu)
    
    if (length(motifs)) {
      for (i in seq_len(nrow(seq_mat))) {
        ch  <- seq_mat[i, ]
        ac  <- ch %in% c("A","C","G","T","a","c","g","t")
        idx <- which(ac)
        if (!length(idx)) next
        
        subj <- Biostrings::DNAString(paste0(toupper(ch[idx]), collapse = ""))
        for (mot in motifs) {
          hits <- Biostrings::matchPattern(mot, subj)
          if (!length(hits)) next
          hs <- Biostrings::start(hits)
          he <- hs + Biostrings::width(hits) - 1L
          for (h in seq_along(hs)) {
            rep_mask[idx[hs[h]:he[h]]] <- TRUE
          }
        }
      }
    }
  }
  
  # 3) combined mask
  anchor_mask    <- hpcol_mask | rep_mask
  cols_to_remove <- which(anchor_mask)
  
  extend_removal <- function(seq_mat, anchor_mask, cols_to_remove) {
    for (col in which(anchor_mask)) {
      # forward
      scan_col <- col + 1L
      while (scan_col <= ncol(seq_mat)) {
        if (scan_col %in% cols_to_remove) {
          scan_col <- scan_col + 1L
          next
        }
        gap_frac <- sum(seq_mat[, scan_col] == "-") / nrow(seq_mat)
        if (gap_frac <= gap_fraction_cutoff) break
        cols_to_remove <- union(cols_to_remove, scan_col)
        scan_col <- scan_col + 1L
      }
      # backward
      scan_col <- col - 1L
      while (scan_col >= 1L) {
        if (scan_col %in% cols_to_remove) {
          scan_col <- scan_col - 1L
          next
        }
        gap_frac <- sum(seq_mat[, scan_col] == "-") / nrow(seq_mat)
        if (gap_frac <= gap_fraction_cutoff) break
        cols_to_remove <- union(cols_to_remove, scan_col)
        scan_col <- scan_col - 1L
      }
    }
    sort(unique(cols_to_remove))
  }
  
  cols_to_remove <- extend_removal(seq_mat, anchor_mask, cols_to_remove)
  
  original_cols  <- ncol(seq_mat)
  filtered_cols  <- length(cols_to_remove)
  remaining_cols <- original_cols - filtered_cols
  
  if (length(cols_to_remove) > 0) {
    filtered_alignment <- seq_mat[, -cols_to_remove, drop = FALSE]
    filtered_phred     <- phred_mat[, -cols_to_remove, drop = FALSE]
  } else {
    filtered_alignment <- seq_mat
    filtered_phred     <- phred_mat
  }
  
  filtered_strings <- apply(filtered_alignment, 1, paste0, collapse = "")
  filtered_quality <- vapply(seq_len(nrow(filtered_phred)), function(i) {
    q <- filtered_phred[i, ]
    q[is.na(q)] <- 0L
    rawToChar(as.raw(q + 33L))
  }, character(1))
  
  stopifnot(all(nchar(filtered_strings) == nchar(filtered_quality)))
  
  fq_header   <- as.character(ShortRead::id(fq))
  fq_filtered <- ShortRead::ShortReadQ(
    sread   = Biostrings::DNAStringSet(filtered_strings),
    quality = Biostrings::BStringSet(filtered_quality),
    id      = Biostrings::BStringSet(fq_header)
  )
  
  output_base <- sub("\\.fastq$", "", basename(fastq_path), ignore.case = TRUE)
  output_file <- file.path(fastq_dir, paste0(output_base, output_suffix, ".fastq"))
  if (file.exists(output_file)) file.remove(output_file)
  ShortRead::writeFastq(fq_filtered, output_file, compress = FALSE)
  
  cbind.data.frame(
    Original        = original_cols,
    Removed         = filtered_cols,
    Remaining       = remaining_cols,
    Percent_Removed = round(100 * filtered_cols / original_cols),
    Output          = output_file
  )
}
