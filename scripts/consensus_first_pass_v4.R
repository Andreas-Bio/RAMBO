generate_consensus_with_export <- function(
    results_cluster_i,
    fastq_dir,
    min_cluster_size = 6,
    phred_threshold = 15,
    iupac_min_count = 3,
    iupac_min_prop = 0.2,
    output_suffix = "_consensus",
    max_pairs = 8^10,
    seed = NULL,
    pairwise_subset_max_inliers = 200L
){
  # --- compile scorer once ---
  if (!exists("hamming_aligned_dual", mode = "function")) {
    Rcpp::cppFunction('
  double hamming_aligned_dual(const std::string& a, const std::string& b, bool ambig) {
    if (a.size() != b.size()) return NA_REAL;
    const size_t n = a.size();
    auto maskAmbig = [](char c)->unsigned {
      switch(std::toupper(c)) {
        case \'A\': return 1u; case \'C\': return 2u; case \'G\': return 4u; case \'T\': return 8u;
        case \'R\': return 1u|4u; case \'Y\': return 2u|8u; case \'S\': return 2u|4u; case \'W\': return 1u|8u;
        case \'K\': return 4u|8u; case \'M\': return 1u|2u; case \'B\': return 2u|4u|8u; case \'D\': return 1u|4u|8u;
        case \'H\': return 1u|2u|8u; case \'V\': return 1u|2u|4u; case \'N\': return 1u|2u|4u|8u;
        default: return 0u;
      }
    };
    auto isACGT = [](char c)->bool {
      switch(std::toupper(c)) { case \'A\': case \'C\': case \'G\': case \'T\': return true; default: return false; }
    };
    unsigned mism = 0u;
    for (size_t i = 0; i < n; ++i) {
      const char p = a[i], q = b[i];
      const bool gp = (p == \'-\'), gq = (q == \'-\');
      if (gp && gq) continue;
      if (gp || gq) { ++mism; continue; }
      if (ambig) {
        unsigned mp = maskAmbig(p), mq = maskAmbig(q);
        if ((mp & mq) == 0u) ++mism;
      } else {
        if (!(isACGT(p) && isACGT(q) && std::toupper(p) == std::toupper(q))) ++mism;
      }
    }
    return 100.0 * static_cast<double>(mism) / static_cast<double>(n);
  }', depends = "Rcpp")
  }
  
  file_stem <- names(results_cluster_i)
  cluster_table_all <- results_cluster_i[[1]]
  fastq_path <- file.path(fastq_dir, paste0(file_stem, ".fastq"))
  
  fq_obj <- ShortRead::readFastq(fastq_path)
  dna_aligned_all_reads <- ShortRead::sread(fq_obj)
  qual_char_by_read <- Biostrings::quality(fq_obj) |> Biostrings::quality() |> as.character()
  read_ids_all <- ShortRead::id(fq_obj) |> as.character()
  
  cluster_table_inliers_only <- dplyr::filter(cluster_table_all, cluster_id != 0)
  if (nrow(cluster_table_inliers_only) == 0) {
    message("All reads labeled as noise, using fallback single cluster.")
    cluster_table_inliers_only <- cluster_table_all
    cluster_table_inliers_only$cluster_id <- 1L
  }
  
  cluster_list_split <- cluster_table_inliers_only |>
    dplyr::group_by(cluster_id) |>
    dplyr::group_split()
  
  summary_rows_list <- vector("list", length(cluster_list_split))
  out_row_idx <- 0L
  
  for (cluster_frame_current in cluster_list_split) {
    cluster_identifier <- unique(cluster_frame_current$cluster_id)
    
    read_index_for_cluster <- match(cluster_frame_current$sequence, read_ids_all) |> stats::na.omit() |> as.integer()
    if (!length(read_index_for_cluster)) next
    
    dna_aligned_cluster_reads <- dna_aligned_all_reads[read_index_for_cluster]
    qual_char_cluster_reads <- qual_char_by_read[read_index_for_cluster]
    read_ids_cluster <- read_ids_all[read_index_for_cluster]
    
    inlier_mask <- !cluster_frame_current$is_outlier
    dna_aligned_inlier_sequences <- dna_aligned_cluster_reads[inlier_mask]
    qual_char_inlier_sequences <- qual_char_cluster_reads[inlier_mask]
    read_ids_inlier <- read_ids_cluster[inlier_mask]
    
    read_ids_outlier <- if (any(!inlier_mask)) paste0(read_ids_cluster[!inlier_mask], "_OUTLIER") else character(0)
    dna_aligned_outlier_sequences_chr <- if (any(!inlier_mask)) as.character(dna_aligned_cluster_reads[!inlier_mask]) else character(0)
    
    if (length(dna_aligned_inlier_sequences) < min_cluster_size) next
    message("Cluster ", cluster_identifier, " final n_seq: ", length(dna_aligned_inlier_sequences))
    
    # consensus
    consensus_sequence_final <- process_consensus_for_cluster(
      dna_aligned_inlier_sequences,
      qual_char_inlier_sequences,
      min_phred       = phred_threshold,
      iupac_min_count = iupac_min_count,
      iupac_min_prop  = iupac_min_prop
    )
    
    # stats
    seqs_chr <- as.character(dna_aligned_inlier_sequences)
    n_inliers <- length(seqs_chr)
    if (!is.null(seed)) set.seed(seed)
    idx_subset <- if (n_inliers > pairwise_subset_max_inliers) sort(sample.int(n_inliers, pairwise_subset_max_inliers, replace = FALSE)) else seq_len(n_inliers)
    seqs_subset_chr <- seqs_chr[idx_subset]
    
    pair_index_matrix <- utils::combn(length(seqs_subset_chr), 2)
    pairwise_distances_percent <- vapply(
      seq_len(ncol(pair_index_matrix)),
      function(k) hamming_aligned_dual(seqs_subset_chr[pair_index_matrix[1, k]],
                                       seqs_subset_chr[pair_index_matrix[2, k]],
                                       FALSE),
      numeric(1)
    )
    imd_percent <- stats::median(pairwise_distances_percent, na.rm = TRUE)
    
    to_consensus_distances_percent <- vapply(
      seq_len(n_inliers),
      function(i) hamming_aligned_dual(seqs_chr[i], consensus_sequence_final, TRUE),
      numeric(1)
    )
    iad_percent <- mean(to_consensus_distances_percent, na.rm = TRUE)
    mdc_percent <- stats::median(to_consensus_distances_percent, na.rm = TRUE)
    
    # per-cluster FASTA (IDs unchanged; only outliers get suffix)
    fasta_per_cluster <- c(
      dna_aligned_inlier_sequences,
      Biostrings::DNAStringSet(dna_aligned_outlier_sequences_chr)
    )
    names(fasta_per_cluster) <- c(read_ids_inlier, read_ids_outlier)
    Biostrings::writeXStringSet(fasta_per_cluster, file.path(fastq_dir, sprintf("%s_CL%s.fasta", file_stem, cluster_identifier)))
    
    out_row_idx <- out_row_idx + 1L
    summary_rows_list[[out_row_idx]] <- list(
      cluster_id                      = cluster_identifier,
      consensus_sequence              = consensus_sequence_final,
      num_sequences                   = length(dna_aligned_inlier_sequences),
      intra_cluster_avg_divergence    = iad_percent,
      intra_cluster_median_divergence = imd_percent,
      mdc_percent                     = mdc_percent,
      num_iupac_codes                 = stringr::str_count(consensus_sequence_final, "[RYKMSWBDHVN]"),
      ungapped_length                 = nchar(gsub("-", "", consensus_sequence_final))
    )
  }
  
  if (out_row_idx == 0L) return(NULL)
  summary_rows_list <- summary_rows_list[seq_len(out_row_idx)]
  
  summary_dataframe <- summary_rows_list |>
    lapply(as.data.frame) |>
    do.call(rbind, args = _)
  
  summary_dataframe$file <- file_stem
  summary_dataframe$cluster_name <- paste0(
    "CL", summary_dataframe$cluster_id,
    "_IAD", round(summary_dataframe$intra_cluster_avg_divergence, 4),
    "_IMD", round(summary_dataframe$intra_cluster_median_divergence, 2),
    "_NIUPAC", summary_dataframe$num_iupac_codes,
    "_NSEQ", summary_dataframe$num_sequences,
    "_LEN", summary_dataframe$ungapped_length
  )
  
  consensus_fasta_all <- Biostrings::DNAStringSet(summary_dataframe$consensus_sequence) |>
    stats::setNames(summary_dataframe$cluster_name)
  
  out_fasta_path <- file.path(fastq_dir, paste0(file_stem, output_suffix, ".fasta"))
  Biostrings::writeXStringSet(consensus_fasta_all, out_fasta_path)
  message("Wrote consensus FASTA to: ", out_fasta_path)
  
  summary_dataframe
}


process_consensus_for_cluster <- function(
    cluster_seqs,
    cluster_qual_char_lists,
    min_phred,
    iupac_min_count,
    iupac_min_prop,
    # Homopolymer-Parameter analog zu filter_homopolymer_alignment
    hp_min_len             = 6L,
    hp_gap_fraction_cutoff = 0.2,
    hp_homopolymer_threshold = 0.5
){
  # ggf. C++ Homopolymer-Scanner nachladen
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
              run_base = "";
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
  
  # Alignment in Matrixform
  seq_matrix <- do.call(rbind, strsplit(as.character(cluster_seqs), ""))
  n_reads    <- nrow(seq_matrix)
  seq_len    <- ncol(seq_matrix)
  
  # Homopolymer-Matrix (gleich wie im Filter)
  hpmat <- get_hpmatrix_cpp(seq_matrix, min_len = hp_min_len)
  
  # Basismaske: Spalten mit vielen HP-Hits
  hpcol_mask_core <- colSums(hpmat) / n_reads >= hp_homopolymer_threshold
  hp_cols_core    <- which(hpcol_mask_core)
  
  # gleiche Erweiterungslogik wie in filter_homopolymer_alignment
  extend_removal <- function(seq_mat, hpcol_mask, cols_to_remove,
                             gap_fraction_cutoff = 0.2){
    for (col in which(hpcol_mask)) {
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
  
  if (length(hp_cols_core)) {
    hp_cols_all <- extend_removal(
      seq_mat           = seq_matrix,
      hpcol_mask        = hpcol_mask_core,
      cols_to_remove    = hp_cols_core,
      gap_fraction_cutoff = hp_gap_fraction_cutoff
    )
  } else {
    hp_cols_all <- integer(0L)
  }
  hp_mask <- logical(seq_len)
  hp_mask[hp_cols_all] <- TRUE
  
  # Qualitäten nach PHRED
  qual_int_lists <- lapply(cluster_qual_char_lists, function(chars) {
    utf8ToInt(paste(chars, collapse = "")) - 33
  })
  qual_matrix <- do.call(rbind, qual_int_lists)
  
  consensus <- character(seq_len)
  
  for (i in seq_len(seq_len)) {
    bases <- seq_matrix[, i]
    quals <- qual_matrix[, i]
    
    is_gap <- bases == "-"
    is_nuc <- bases %in% c("A", "C", "G", "T")
    
    # PHRED-Filter nur für Nukleotide, Gaps immer behalten
    valid_idx <- which((quals >= min_phred & is_nuc) | is_gap)
    
    # Fallback bei sehr wenig Signal
    if (length(valid_idx) < 3) {
      ord <- order(quals, decreasing = TRUE)
      valid_idx <- ord[seq_len(min(3L, length(ord)))]
    }
    
    filtered_bases <- bases[valid_idx]
    base_counts    <- table(filtered_bases)
    
    nuc_counts <- base_counts[names(base_counts) %in% c("A","C","G","T")]
    if (is.null(nuc_counts)) nuc_counts <- integer(0)
    gap_count <- base_counts["-"]
    gap_total <- ifelse(is.na(gap_count), 0L, as.integer(gap_count))
    nuc_total <- sum(nuc_counts, na.rm = TRUE)
    
    # Homopolymer-Spalte: immer stricter Mehrheits-Call
    if (hp_mask[i]) {
      # nur A/C/G/T/- betrachten
      hp_bases <- intersect(names(base_counts), c("A","C","G","T","-"))
      if (!length(hp_bases)) {
        consensus[i] <- "N"
      } else {
        hp_counts <- base_counts[hp_bases]
        # strikte Mehrheit (bei Gleichstand: erste mit max, wie table-Order)
        winner <- hp_bases[which.max(hp_counts)]
        consensus[i] <- winner
      }
    } else {
      # Standardfall: Gap vs Nukleotide + IUPAC
      if (gap_total > nuc_total) {
        consensus[i] <- "-"
      } else {
        consensus[i] <- resolve_iupac(nuc_counts, iupac_min_count, iupac_min_prop)
      }
    }
  }
  
  paste0(consensus, collapse = "")
}



# Resolve IUPAC based on frequency and count thresholds
resolve_iupac <- function(col_counts, min_count, min_prop) {
  total <- sum(col_counts)
  if (total == 0) return("N")
  
  # Identify bases above thresholds
  above_threshold <- col_counts >= min_count & (col_counts / total) >= min_prop
  bases <- names(col_counts)[above_threshold]
  
  if (length(bases) == 0) {
    # If nothing meets the threshold, return the most frequent base
    bases <- names(col_counts)[which.max(col_counts)]
  }
  
  # Build lookup once, outside this function if reused often
  iupac_keys <- sapply(IUPAC_CODE_MAP, function(b) paste0(sort(strsplit(b, "")[[1]]), collapse = ""))
  unique_keys <- !duplicated(iupac_keys)
  iupac_lookup <- setNames(names(IUPAC_CODE_MAP)[unique_keys], iupac_keys[unique_keys])
  
  # Sort selected bases and generate key
  key <- paste0(sort(bases), collapse = "")
  
  # Map to IUPAC code if possible
  iupac_code <- iupac_lookup[[key]]
  if (is.null(iupac_code)) {
    # Fallback: concatenate bases (not ideal but prevents failure)
    iupac_code <- paste0(bases, collapse = "")
  }
  
  return(iupac_code)
}



