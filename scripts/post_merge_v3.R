merge_and_regenerate_consensus <- function(
    input_fasta_file,
    min_cluster_size        = 6,
    phred_threshold         = 15,
    iupac_min_count         = 3,
    iupac_min_prop          = 0.2,
    max_pairwise_distance   = 0,     # 0 => only identical consensuses merge
    align_suffix            = align_fastq_with_mafft.output_suffix_align,
    iupac_guard             = TRUE,
    iupac_delta_abs         = NULL,  # NULL => use max_pairwise_distance
    iupac_delta_rel         = 0      # e.g. 0.005 allows 0.5% of ungapped length
){
  stopifnot(file.exists(input_fasta_file))
  suppressPackageStartupMessages({
    library(Biostrings); library(ShortRead); library(IRanges); library(stringr)
  })
  
  # --- helpers ---------------------------------------------------------------
  find_aligned_fastq <- function(stem, dir, align_suffix) {
    files <- list.files(dir, pattern = "\\.fastq$", full.names = TRUE)
    base  <- basename(files)
    hits  <- files[startsWith(base, stem) & endsWith(base, paste0(align_suffix, ".fastq"))]
    if (!length(hits)) stop("No aligned FASTQ for stem=", stem,
                            " suffix=", align_suffix, " in ", dir)
    hits[[1]]
  }
  
  # C++ now loaded externally. Compile once with:
  # Rcpp::sourceCpp("hamming_mismatches_aligned_dual.cpp")
  if (!exists("hamming_mismatches_aligned_dual")) {
    Rcpp::cppFunction("
SEXP hamming_mismatches_aligned_dual(Rcpp::CharacterVector seqs,
                                     bool ambig = false,
                                     bool return_percent = false) {
  int n = seqs.size(); if (n == 0) Rcpp::stop(\"Need >=1 seq\");
  std::vector<std::string> s(n);
  size_t L = Rcpp::as<std::string>(seqs[0]).size();
  for (int i = 0; i < n; ++i) {
    s[i] = Rcpp::as<std::string>(seqs[i]);
    if (s[i].size() != L) Rcpp::stop(\"All sequences equal length.\");
  }
  unsigned mask[256]; std::fill(mask, mask+256, 0u);
  auto setmask=[&](char c,unsigned m){
    mask[(unsigned char)std::toupper((unsigned char)c)]=m;
    mask[(unsigned char)std::tolower((unsigned char)c)]=m;
  };
  setmask('A',1u); setmask('C',2u); setmask('G',4u); setmask('T',8u);
  setmask('R',1u|4u); setmask('Y',2u|8u); setmask('S',2u|4u); setmask('W',1u|8u);
  setmask('K',4u|8u); setmask('M',1u|2u); setmask('B',2u|4u|8u); setmask('D',1u|4u|8u);
  setmask('H',1u|2u|8u); setmask('V',1u|2u|4u); setmask('N',1u|2u|4u|8u);
  auto is_acgt=[](unsigned char u)->bool { return (u=='A'||u=='C'||u=='G'||u=='T'); };

  Rcpp::NumericMatrix D(n,n); for (int i=0;i<n;++i) D(i,i)=0.0;
  for (int i=0;i<n-1;++i){
    const std::string &a=s[i];
    for (int j=i+1;j<n;++j){
      const std::string &b=s[j];
      size_t mis=0;
      for (size_t k=0;k<L;++k){
        char ca=a[k], cb=b[k];
        if (ca=='-' && cb=='-') continue;
        if (ca=='-' || cb=='-') { ++mis; continue; }
        unsigned char ua=std::toupper((unsigned char)ca),
                      ub=std::toupper((unsigned char)cb);
        if (ambig) { if ((mask[ua] & mask[ub])==0u) ++mis; }
        else       { if (!(is_acgt(ua)&&is_acgt(ub)&&ua==ub)) ++mis; }
      }
      double val = return_percent ? (L?(100.0*(double)mis/(double)L):NA_REAL) : (double)mis;
      D(i,j)=D(j,i)=val;
    }
  }
  if (n==2) return return_percent?Rcpp::wrap(D(0,1)):Rcpp::wrap((int)D(0,1));
  return D;
}", plugins = "cpp11")
  }
  
  .as_square <- function(Draw, nm) if (is.matrix(Draw)) { dimnames(Draw) <- list(nm,nm); Draw } else {
    matrix(c(0,Draw,Draw,0), 2, 2, dimnames = list(nm,nm))
  }
  
  .iupac_mask <- c(A=1L,C=2L,G=4L,T=8L,R=5L,Y=10L,S=6L,W=9L,K=12L,M=3L,B=14L,D=13L,H=11L,V=7L,N=15L,"-"=0L)
  homopoly_mask_iupac_allow_gaps <- function(s, min_len = 6L) {
    x <- strsplit(toupper(s), "", fixed = TRUE)[[1]]; n <- length(x); m <- rep(FALSE, n)
    i <- 1L
    while (i <= n) {
      ai <- x[i]; if (!(ai %in% c("A","C","G","T"))) { i <- i + 1L; next }
      anchor <- .iupac_mask[ai]; j <- i + 1L; run <- 1L
      while (j <= n) {
        aj <- x[j]; if (aj == "-") { j <- j + 1L; next }
        if (bitwAnd(anchor, .iupac_mask[aj]) == 0L) break
        run <- run + 1L; j <- j + 1L
      }
      if (run >= min_len) m[i:(j-1L)] <- TRUE
      i <- j
    }
    m
  }
  
  trim_and_dist_for_merge <- function(cons_set, hp_min = 6L, di_rep = 4L, tri_rep = 4L) {
    L <- unique(width(cons_set)); stopifnot(length(L)==1L); L <- as.integer(L)
    del <- rep(FALSE, L)
    
    for (s in as.character(cons_set)) del[ homopoly_mask_iupac_allow_gaps(s, hp_min) ] <- TRUE
    
    bases <- c("A","C","G","T")
    di  <- as.vector(outer(bases, bases, paste0))
    tri <- as.vector(outer(outer(bases, bases, paste0), bases, paste0))
    dinu <- vapply(di,  function(m) paste0(rep(m,  di_rep),  collapse=""), "")
    triu <- vapply(tri, function(m) paste0(rep(m, tri_rep), collapse=""), "")
    
    for (s in as.character(cons_set)) {
      ch <- strsplit(s, "", fixed=TRUE)[[1]]
      ac <- ch %in% c("A","C","G","T","a","c","g","t")
      idx <- which(ac); if (!length(idx)) next
      subj <- DNAString(paste0(toupper(ch[idx]), collapse=""))
      for (mot in c(dinu, triu)) {
        hits <- matchPattern(mot, subj); if (!length(hits)) next
        hs <- start(hits); he <- hs + width(hits) - 1L
        for (h in seq_along(hs)) del[idx[hs[h]]:idx[he[h]]] <- TRUE
      }
    }
    
    keep <- !del
    trimmed <- vapply(as.character(cons_set), function(s) {
      ch <- strsplit(s, "", fixed=TRUE)[[1]]
      if (!any(keep)) "" else paste0(ch[keep], collapse="")
    }, "", USE.NAMES = FALSE) |> DNAStringSet() |> `names<-`(names(cons_set))
    
    nm <- names(trimmed)
    Draw <- hamming_mismatches_aligned_dual(as.character(trimmed), ambig=TRUE, return_percent=FALSE)
    list(trimmed = trimmed, D_counts = .as_square(Draw, nm), kept_length = sum(keep))
  }
  
  process_consensus_for_cluster <- function(cluster_seqs, cluster_qual_char_lists,
                                            min_phred, iupac_min_count, iupac_min_prop) {
    seq_mat  <- do.call(rbind, strsplit(as.character(cluster_seqs), ""))
    qual_mat <- do.call(rbind, lapply(cluster_qual_char_lists, \(q) utf8ToInt(paste(q, collapse="")) - 33))
    L <- ncol(seq_mat); cons <- character(L)
    for (i in seq_len(L)) {
      b <- seq_mat[,i]; q <- qual_mat[,i]
      is_gap <- b=="-"; is_nuc <- b %in% c("A","C","G","T")
      valid <- which((q >= min_phred & is_nuc) | is_gap)
      if (length(valid) < 3) valid <- order(q, decreasing=TRUE)[seq_len(min(3, length(q)))]
      tbl <- table(b[valid]); nuc <- tbl[names(tbl)%in%c("A","C","G","T")]
      gap <- ifelse(is.na(tbl["-"]), 0, tbl["-"]); nuc_n <- sum(nuc, na.rm=TRUE)
      if (gap > nuc_n) { cons[i] <- "-" } else {
        total <- sum(nuc); if (total==0) { cons[i] <- "N"; next }
        keep <- nuc >= iupac_min_count & (nuc/total) >= iupac_min_prop
        bases <- names(nuc)[keep]; if (!length(bases)) bases <- names(nuc)[which.max(nuc)]
        keys <- sapply(IUPAC_CODE_MAP, \(b) paste0(sort(strsplit(b, "")[[1]]), collapse=""))
        lk   <- setNames(names(IUPAC_CODE_MAP)[!duplicated(keys)], keys[!duplicated(keys)])
        cons[i] <- lk[[ paste0(sort(bases), collapse="") ]] %||% paste0(bases, collapse="")
      }
    }
    paste0(cons, collapse="")
  }
  
  `%||%` <- function(a,b) if (is.null(a)) b else a
  strip_outlier <- function(x) sub("_OUTLIER$", "", x, perl = TRUE)
  .count_iupac <- function(s) stringr::str_count(s, "[RYKMSWBDHVN]")
  
  pairwise_IMD_percent <- function(inlier_seq_char, cap = 200L) {
    n <- length(inlier_seq_char); if (n < 2) return(NA_real_)
    idx <- if (n > cap) sort(sample.int(n, cap)) else seq_len(n)
    ss  <- inlier_seq_char[idx]; L <- unique(nchar(ss)); stopifnot(length(L)==1L)
    D   <- hamming_mismatches_aligned_dual(ss, ambig=FALSE, return_percent=FALSE)
    stats::median((if (is.matrix(D)) D[upper.tri(D)] else as.numeric(D)) / L * 100, na.rm=TRUE)
  }
  
  # --- I/O -------------------------------------------------------------------
  stem <- sub("_final$", "", tools::file_path_sans_ext(basename(input_fasta_file)))
  dir_out <- dirname(input_fasta_file)
  
  # clean old merged outputs
  allf <- list.files(dir_out, full.names = TRUE)
  bn <- basename(allf)
  # stringr-basiertes Escaping ohne gsub; Backslash als ASCII 92
  escape_regex_fixed <- function(s) {
    specials <- c(intToUtf8(92), "[", "]", "{", "}", "(", ")", "+", "*", "^", "$", ".", "|", "?", "-")
    out <- s
    for (ch in specials) {
      out <- stringr::str_replace_all(out, stringr::fixed(ch), paste0(intToUtf8(92), ch))
    }
    out
  }
  rx <- escape_regex_fixed(stem)
  dead <- stringr::str_detect(bn, stringr::regex(paste0("^", rx, "_CL[0-9_]+_merged[.]fasta$"))) |
    stringr::str_detect(bn, stringr::regex(paste0("^", rx, "_final_merged[.]fasta$")))
  if (any(dead)) suppressWarnings(unlink(allf[dead], force = TRUE))
  
  per_cluster_fasta <- list.files(dir_out,
                                  pattern = paste0("^", rx, "_CL[0-9]+[.]fasta$"), full.names = TRUE)
  if (!length(per_cluster_fasta)) stop("No per-cluster FASTA in ", dir_out)
  
  cons_all <- readDNAStringSet(input_fasta_file)
  if (!length(cons_all)) stop("Empty consensus file: ", input_fasta_file)
  
  if (length(cons_all) == 1L) {
    out1 <- file.path(dir_out, paste0(stem, "_final_merged.fasta"))
    writeXStringSet(cons_all, out1)
    nm <- names(cons_all)
    return(list(cluster_names = nm,
                membership_map = list(setNames(nm, nm)),
                initial_consensus_dist_matrix = matrix(0,1,1, dimnames=list(nm,nm))))
  }
  
  cn   <- names(cons_all)
  cid  <- stringr::str_extract(cn, "^CL[0-9]+")
  nseq <- as.integer(stringr::str_extract(cn, "(?<=_NSEQ)[0-9]+"))
  keep <- !is.na(cid) & nseq >= min_cluster_size
  if (!any(keep)) stop("No consensus meets min_cluster_size=", min_cluster_size)
  
  cons_tbl <- data.frame(
    consensus_name = cn[keep],
    consensus_seq  = as.character(cons_all)[keep],
    cluster_id     = cid[keep],
    nseq           = nseq[keep],
    stringsAsFactors = FALSE
  )
  
  # map cluster id -> fasta path, pre-read to avoid repeated I/O
  pc_ids <- stringr::str_extract(basename(per_cluster_fasta), "CL[0-9]+")
  mm     <- match(pc_ids, cons_tbl$cluster_id)
  per_cluster_fasta <- per_cluster_fasta[!is.na(mm)]
  pc_ids            <- pc_ids[!is.na(mm)]
  
  # cache reads per cluster
  reads_by_id     <- setNames(vector("list", length(pc_ids)), pc_ids)
  reads_in_by_id  <- setNames(vector("list", length(pc_ids)), pc_ids)
  names_core_by_id<- setNames(vector("list", length(pc_ids)), pc_ids)
  for (i in seq_along(pc_ids)) {
    rd <- readDNAStringSet(per_cluster_fasta[i])
    reads_by_id[[i]] <- rd
    in_mask <- !endsWith(names(rd), "_OUTLIER")
    rd_in   <- rd[in_mask]
    reads_in_by_id[[i]] <- rd_in
    names_core_by_id[[i]] <- sub("_OUTLIER$", "", names(rd_in), perl=TRUE)
  }
  
  id_to_idx <- setNames(seq_along(pc_ids), pc_ids)
  
  aligned_fastq_path <- find_aligned_fastq(stem, dir_out, align_suffix)
  fq <- ShortRead::readFastq(aligned_fastq_path)
  ali_ids   <- as.character(ShortRead::id(fq))
  ali_qchar <- Biostrings::quality(fq) |> Biostrings::quality() |> as.character()
  ali_index <- setNames(seq_along(ali_ids), ali_ids)
  
  # --- merge groups by trimmed-count distances -------------------------------
  cons_set <- DNAStringSet(cons_tbl$consensus_seq) |> `names<-`(cons_tbl$cluster_id)
  td       <- trim_and_dist_for_merge(cons_set, hp_min = 6L, di_rep = 4L, tri_rep = 4L)
  Dcounts  <- td$D_counts
  
  groups <- if (max_pairwise_distance < 0 || length(cons_set) < 2L) {
    setNames(as.list(names(cons_set)), names(cons_set))
  } else {
    hc   <- hclust(as.dist(Dcounts), method = "single")
    lab  <- cutree(hc, h = max_pairwise_distance)
    sp   <- split(names(cons_set), lab)
    lapply(sp, sort) |> `names<-`(vapply(sp, `[[`, character(1), 1))
  }
  
  .build_merge_with_iupac_guard <- function(
    merge_ids, min_cluster_size,
    ali_index, ali_qchar,
    phred_threshold, iupac_min_count, iupac_min_prop,
    iupac_delta_abs, iupac_delta_rel
  ){
    acc_ids <- character(0)
    pooled_reads <- DNAStringSet()
    pooled_in    <- DNAStringSet()
    current_cons <- NULL
    current_iup  <- 0L
    current_ulen <- 0L
    
    for (id in merge_ids) {
      j  <- id_to_idx[[id]]
      rd <- reads_by_id[[j]]
      rd_in <- reads_in_by_id[[j]]
      if (length(rd_in) < min_cluster_size && length(acc_ids) == 0L) next
      
      rd_in_new <- c(pooled_in, rd_in)
      if (length(rd_in_new) < min_cluster_size) next
      
      rn_core <- c(unlist(names_core_by_id[match(acc_ids, pc_ids)]), names_core_by_id[[j]])
      idx <- unname(ali_index[rn_core])
      if (anyNA(idx)) next
      
      cons_new <- process_consensus_for_cluster(
        cluster_seqs            = rd_in_new,
        cluster_qual_char_lists = ali_qchar[idx],
        min_phred       = phred_threshold,
        iupac_min_count = iupac_min_count,
        iupac_min_prop  = iupac_min_prop
      )
      n_iup_new <- .count_iupac(cons_new)
      ulen_new  <- nchar(gsub("-", "", cons_new))
      
      abs_limit <- if (is.null(iupac_delta_abs)) as.integer(max(0, max_pairwise_distance)) else as.integer(iupac_delta_abs)
      rel_limit <- as.integer(ceiling(iupac_delta_rel * max(1L, ulen_new)))
      allow_inc <- max(abs_limit, rel_limit)
      
      inc <- n_iup_new - current_iup
      if (length(acc_ids) == 0L) {
        acc_ids      <- c(acc_ids, id)
        pooled_reads <- rd
        pooled_in    <- rd_in
        current_cons <- cons_new
        current_iup  <- n_iup_new
        current_ulen <- ulen_new
      } else if (inc <= allow_inc) {
        acc_ids      <- c(acc_ids, id)
        pooled_reads <- c(pooled_reads, rd)
        pooled_in    <- rd_in_new
        current_cons <- cons_new
        current_iup  <- n_iup_new
        current_ulen <- ulen_new
      } else {
        # reject this candidate
        next
      }
    }
    list(accepted_ids = acc_ids,
         pooled_reads = pooled_reads,
         consensus    = current_cons,
         n_iupac      = current_iup,
         ulen         = current_ulen)
  }
  
  merged_cons  <- list()
  merged_mems  <- list()
  skipped_small <- 0L
  
  for (seed_id in names(groups)) {
    merge_ids <- groups[[seed_id]]
    ord <- order(Dcounts[seed_id, merge_ids], na.last = NA)
    merge_ids <- merge_ids[ord]
    
    built <- if (!iupac_guard) {
      member_idx  <- unname(id_to_idx[merge_ids])
      pooled_reads<- do.call(c, reads_by_id[member_idx])
      list(
        accepted_ids = merge_ids,
        pooled_reads = pooled_reads,
        consensus    = NULL,
        n_iupac      = NA_integer_,
        ulen         = NA_integer_
      )
    } else {
      .build_merge_with_iupac_guard(
        merge_ids           = merge_ids,
        min_cluster_size    = min_cluster_size,
        ali_index           = ali_index,
        ali_qchar           = ali_qchar,
        phred_threshold     = phred_threshold,
        iupac_min_count     = iupac_min_count,
        iupac_min_prop      = iupac_min_prop,
        iupac_delta_abs     = iupac_delta_abs,
        iupac_delta_rel     = iupac_delta_rel
      )
    }
    
    if (!length(built$accepted_ids)) { skipped_small <- skipped_small + 1L; next }
    
    # compute consensus if not already available
    pooled_reads <- built$pooled_reads
    in_mask <- !endsWith(names(pooled_reads), "_OUTLIER")
    pooled_in <- pooled_reads[in_mask]
    rn_core <- sub("_OUTLIER$", "", names(pooled_in), perl=TRUE)
    idx <- unname(ali_index[rn_core])
    
    merged_seq <- if (is.null(built$consensus)) {
      process_consensus_for_cluster(
        cluster_seqs            = pooled_in,
        cluster_qual_char_lists = ali_qchar[idx],
        min_phred       = phred_threshold,
        iupac_min_count = iupac_min_count,
        iupac_min_prop  = iupac_min_prop
      )
    } else built$consensus
    
    in_chr <- as.character(pooled_in)
    Lr <- unique(nchar(in_chr)); if (length(Lr) != 1L) stop("Inlier reads must share the same alignment length.")
    imd <- pairwise_IMD_percent(in_chr, cap = 200L)
    d_to_cons <- vapply(seq_along(in_chr), function(i)
      hamming_mismatches_aligned_dual(c(in_chr[i], merged_seq), ambig=TRUE, return_percent=FALSE),
      numeric(1))
    mdc <- stats::median(100 * d_to_cons / Lr, na.rm = TRUE)
    n_iup <- stringr::str_count(merged_seq, "[RYKMSWBDHVN]")
    ulen  <- nchar(gsub("-", "", merged_seq))
    nin   <- length(in_chr)
    
    core <- paste0("CL", paste0(sub("^CL", "", sort(built$accepted_ids)), collapse = "_"))
    label <- paste0(core,
                    "_IMD", formatC(imd, digits=2, format="f"),
                    "_MDC", formatC(mdc, digits=2, format="f"),
                    "_NIUPAC", n_iup,
                    "_NSEQ", nin,
                    "_LEN", ulen)
    
    out_reads <- file.path(dir_out, paste0(stem, "_", core, "_merged.fasta"))
    Biostrings::writeXStringSet(pooled_reads, out_reads)
    
    merged_cons[[label]] <- merged_seq
    merged_mems[[label]] <- built$accepted_ids
  }
  
  out_fasta <- if (length(merged_cons)) {
    DNAStringSet(unlist(merged_cons, use.names = FALSE)) |> `names<-`(names(merged_cons))
  } else {
    DNAStringSet(cons_tbl$consensus_seq) |> `names<-`(cons_tbl$cluster_id)
  }
  
  final_path <- file.path(dir_out, paste0(stem, "_final_merged.fasta"))
  if (file.exists(final_path)) suppressWarnings(unlink(final_path, force = TRUE))
  writeXStringSet(out_fasta, final_path)
  
  message("Sample ", stem, ": ",
          nrow(cons_tbl), " input clusters, ",
          length(merged_cons), " merged groups written, ",
          skipped_small, " skipped for size.")
  
  list(
    cluster_names = names(merged_cons),
    membership_map = merged_mems,
    initial_consensus_dist_matrix = Dcounts
  )
}


collect_final_merged_consensus <- function(
    input_dir      = ".",
    suffix_pattern = "_final_merged\\.fasta$",
    input_suffix   = NULL
) {
  if (is.null(input_suffix)) {
    if (!exists("input_suffix", envir = .GlobalEnv, inherits = TRUE)) {
      stop("input_suffix must be provided or exist in the global environment")
    }
    input_suffix <- get("input_suffix", envir = .GlobalEnv, inherits = TRUE)
  }
  
  fasta_files <- list.files(input_dir, pattern = suffix_pattern, full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    stop("No *_final_merged.fasta files found in: ", input_dir)
  }
  
  all_metadata <- lapply(fasta_files, function(f) {
    seqs <- Biostrings::readDNAStringSet(f)
    if (length(seqs) == 0) return(NULL)
    
    file_base <- basename(f)
    
    sample_name <- sub(paste0("_", input_suffix, ".*"), "", file_base)
    
    if (identical(sample_name, file_base)) {
      warning("Could not extract sample name from file: ", file_base)
      return(NULL)
    }
    
    data.frame(
      sample     = sample_name,
      cluster_id = stringr::str_extract(names(seqs), "^CL\\d+(?:_\\d+)*"),
      IMD        = as.numeric(stringr::str_extract(names(seqs), "(?<=_IMD)[0-9.]+")),
      MDC        = as.numeric(stringr::str_extract(names(seqs), "(?<=_MDC)[0-9.]+")),
      NIUPAC     = as.integer(stringr::str_extract(names(seqs), "(?<=_NIUPAC)\\d+")),
      NSEQ       = as.integer(stringr::str_extract(names(seqs), "(?<=_NSEQ)\\d+")),
      LEN        = as.integer(stringr::str_extract(names(seqs), "(?<=_LEN)\\d+")),
      stringsAsFactors = FALSE
    )
  })
  
  all_metadata_df <- dplyr::bind_rows(all_metadata)
  return(all_metadata_df)
}



