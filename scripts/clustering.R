detect_rare_minor_dual <- function(
    fq,
    window = detect_rare_minor_dual.window,
    max_gap_frac = detect_rare_minor_dual.max_gap_frac,
    min_frac = detect_rare_minor_dual.min_frac,
    min_reads = detect_rare_minor_dual.min_reads,
    fdr = detect_rare_minor_dual.fdr,
    use_median_bg = detect_rare_minor_dual.use_median_bg
){
  # input
  dna <- Biostrings::DNAStringSet(ShortRead::sread(fq))
  stopifnot(length(unique(Biostrings::width(dna))) == 1L)
  st <- compute_col_stats(dna)
  
  # global column filter
  keep <- which(st$gap_frac <= max_gap_frac & st$depth_acgt > 0)
  if (!length(keep)){
    return(list(minor_table = data.frame(), features = NULL, weights = numeric(0)))
  }
  
  # local background of nonmajor per kept column
  nonmajor_bg <- if (use_median_bg) sliding_stat(st$nonmajor[keep], window, median)
  else               sliding_stat(st$nonmajor[keep], window, mean)
  
  # per column rare minor test with BH
  acgt <- c("A","C","G","T")
  bg_floor <- 5e-4
  rows <- list()
  for (t in seq_along(keep)){
    j <- keep[t]
    n <- st$depth_acgt[j]; if (n <= 0) next
    bg <- max(0, min(0.25, nonmajor_bg[t]))
    k_vec <- as.numeric(st$counts_acgt[, j]); names(k_vec) <- acgt
    maj <- which.max(st$freq_acgt[, j])
    thr_reads <- max(min_reads, ceiling(min_frac * n))
    for (b in setdiff(acgt, acgt[maj])){
      k <- k_vec[b]; frac <- k / n; delta <- frac - pmax(bg, bg_floor)
      if (k < thr_reads || delta < min_frac) next
      p_bg <- max(bg, bg_floor)
      pval <- stats::pbinom(q = k-1, size = n, prob = p_bg, lower.tail = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(pos = j, base = b, k = k, n = n,
                                              frac = frac, bg = bg, delta = delta, p = pval,
                                              stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)){
    return(list(minor_table = data.frame(), features = NULL, weights = numeric(0)))
  }
  
  tab <- do.call(rbind, rows)
  tab$q <- stats::p.adjust(tab$p, method = "BH")
  ord <- order(tab$pos, -tab$delta, tab$q, -tab$k)
  tab <- tab[ord, , drop = FALSE]
  
  keep_idx <- unlist(tapply(seq_len(nrow(tab)), tab$pos, function(ii){
    ii[ tab$q[ii] <= fdr & tab$delta[ii] >= min_frac ]
  }))
  if (!length(keep_idx)){
    return(list(minor_table = data.frame(), features = NULL, weights = numeric(0)))
  }
  
  tab2 <- tab[keep_idx, , drop = FALSE]
  gap_map <- st$gap_frac[tab2$pos]
  
  # base weights before multi base pruning
  w0 <- pmax(0, tab2$delta) * log1p(tab2$n) * (1 - gap_map)
  
  # resolve multiple bases per position in a data driven way
  pos_multi <- as.integer(names(which(table(tab2$pos) > 1L)))
  keep_cols <- rep(TRUE, nrow(tab2))
  if (length(pos_multi)){
    for (p in pos_multi){
      jj <- which(tab2$pos == p)
      if (length(jj) < 2L) next
      B <- sapply(jj, function(j){
        as.integer(as.character(Biostrings::subseq(dna, start = tab2$pos[j], width = 1L)) == tab2$base[j])
      })
      A  <- crossprod(B)
      n1 <- colSums(B)
      U  <- outer(n1, n1, `+`) - A
      Jacc <- suppressWarnings(A / pmax(1L, U)); diag(Jacc) <- 0
      if (any(Jacc > 0)){
        ordj <- order(w0[jj], decreasing = TRUE)
        drop <- jj[ ordj[2:length(jj)] ]
        keep_cols[drop] <- FALSE
        next
      }
      C <- suppressWarnings(stats::cor(B)); C[!is.finite(C)] <- 0
      neg_pairs <- sum(C[lower.tri(C)] < 0)
      K <- min(length(jj), 1 + neg_pairs)
      if (K < length(jj)){
        ordj <- order(w0[jj], decreasing = TRUE)
        drop <- jj[ ordj[(K+1):length(jj)] ]
        keep_cols[drop] <- FALSE
      }
    }
    tab2    <- tab2[keep_cols, , drop = FALSE]
    gap_map <- gap_map[keep_cols]
    w0      <- w0[keep_cols]
  }
  
  # guard if multi base pruning removed everything
  if (!nrow(tab2)){
    return(list(minor_table = data.frame(), features = NULL, weights = numeric(0)))
  }
  
  #purity sharpness idf weighting that is safe for IUPAC
  ids <- as.character(ShortRead::id(fq))
  N   <- length(ids)
  f_pos <- st$freq_acgt[, tab2$pos, drop = FALSE]
  if (!ncol(f_pos)){
    return(list(minor_table = data.frame(), features = NULL, weights = numeric(0)))
  }
  p_sorted <- apply(f_pos, 2, function(v) sort(as.numeric(v), decreasing = TRUE))
  if (!is.matrix(p_sorted)) p_sorted <- matrix(p_sorted, nrow = 4)
  p1 <- p_sorted[1, ]; p2 <- p_sorted[2, ]; p3 <- p_sorted[3, ]
  sep_pos <- p1 - p2
  g23_pos <- pmax(0, p2 - p3)
  g23n    <- ifelse(p2 > 0, g23_pos / p2, 1)
  idf     <- log1p(N / pmax(1, tab2$k))
  
  tab2$weight <- w0 * pmax(1e-6, sep_pos) * (0.5 + 0.5*g23n) * idf
  
  #purity + biallelic-balance + idf weighting, relativ zur globalen Gewichtsverteilung
  
  # build binary feature matrix Z
  J <- nrow(tab2)
  if (J == 0L){
    return(list(minor_table = data.frame(), features = NULL, weights = numeric(0)))
  }
  Z <- matrix(0L, nrow = N, ncol = J)
  colnames(Z) <- paste0("pos", tab2$pos, ":", tab2$base)
  rownames(Z) <- ids
  for (jj in seq_len(J)){
    p <- tab2$pos[jj]; b <- tab2$base[jj]
    col_b <- as.character(Biostrings::subseq(dna, start = p, width = 1L))
    Z[, jj] <- as.integer(col_b == b)
  }
  
  list(minor_table = tab2, features = Z, weights = tab2$weight)
}

compute_col_stats <- function(dna){
  cm <- consensus_matrix_ACGTN_gap(dna)
  counts_acgt <- cm[c("A","C","G","T"), , drop = FALSE]
  depth_acgt  <- colSums(counts_acgt)
  gap_frac    <- as.numeric(cm["-", ] / pmax(1L, colSums(cm)))
  freq_acgt   <- sweep(counts_acgt, 2, pmax(1L, depth_acgt), "/")
  nonmajor    <- 1 - apply(freq_acgt, 2, max)
  list(counts_acgt = counts_acgt, depth_acgt = depth_acgt,
       gap_frac = gap_frac, freq_acgt = freq_acgt, nonmajor = nonmajor)
}

consensus_matrix_ACGTN_gap <- function(dna){
  cm <- Biostrings::consensusMatrix(dna, baseOnly = FALSE, as.prob = FALSE)
  need <- c("A","C","G","T","N","-")
  miss <- setdiff(need, rownames(cm))
  if (length(miss)){
    add <- matrix(0L, nrow = length(miss), ncol = ncol(cm))
    rownames(add) <- miss
    colnames(add) <- colnames(cm)
    cm <- rbind(cm, add)
  }
  cm
}

sliding_stat <- function(x, window = sliding_stat.window, fun = sliding_stat.fun){
  n <- length(x); y <- numeric(n)
  for (i in seq_len(n)){
    a <- max(1L, i-window); b <- min(n, i+window)
    y[i] <- suppressWarnings(fun(x[a:b], na.rm = TRUE))
    if (!is.finite(y[i])) y[i] <- 0
  }
  y
}

hierarchical_column_weights <- function(
    aligned_seqs,
    feature_pos = NULL,              # 1-basierte Positionsvektoren (z B minor_table$pos)
    min_minor_count = 3,
    min_minor_frac  = 0.1,
    min_node_size   = 3,
    n_steps         = 50,
    apply_global_floor = TRUE,
    include_gaps    = TRUE,          # nur für Gap-Filter, Kontrast nur auf ACGT
    step_mode       = "quantile_merge_heights",
    min_branch_frac = 0.1,
    min_branch_size = 10,
    quantile_skew   = 1,
    max_gap_fraction = 0.5,
    min_step_height_frac = 0.01,
    extra_col_frac  = 0.10,          # 10 Prozent Random-Spalten für Baum
    feature_min_for_subsample = 50L  # nur Zusatzspalten wenn Features < 50
){
  n <- length(aligned_seqs)
  if (n == 0L) stop("aligned_seqs is empty")
  
  aln_mat <- do.call(rbind, strsplit(as.character(aligned_seqs), ""))
  aln_len <- ncol(aln_mat)
  
  base_levels <- c("A","C","G","T","N","-")
  base_lookup <- setNames(seq_along(base_levels), base_levels)
  
  # map bases to integer bins
  base_factor_mat <- matrix(NA_integer_, nrow = n, ncol = aln_len)
  for (j in seq_len(aln_len)) {
    base_factor_mat[, j] <- base_lookup[aln_mat[, j]]
  }
  
  # globale Gap-Filter, um extrem gap-lastige Spalten abzuschneiden
  cm_global <- {
    counts <- matrix(0L, nrow = length(base_levels), ncol = aln_len)
    rownames(counts) <- base_levels
    for (j in seq_len(aln_len)) {
      v <- base_factor_mat[, j]
      v <- v[!is.na(v)]
      if (length(v)) {
        tab <- tabulate(v, nbins = length(base_levels))
        counts[, j] <- tab
      }
    }
    counts
  }
  gap_counts <- cm_global["-", ]
  depth_all  <- colSums(cm_global[c("A","C","G","T"), , drop = FALSE])
  gap_frac   <- ifelse(depth_all + gap_counts > 0,
                       gap_counts / (depth_all + gap_counts),
                       0)
  
  keep_cols_global <- which(is.finite(gap_frac) & gap_frac <= max_gap_fraction)
  if (!length(keep_cols_global)) {
    warning("No columns pass global gap filter")
    col_weights <- rep(0, aln_len)
    names(col_weights) <- seq_len(aln_len)
    attr(col_weights, "hit_counts")        <- integer(aln_len)
    attr(col_weights, "avg_cols_per_step") <- rep(NA_real_, aln_len)
    attr(col_weights, "step_heights")      <- NA_real_
    return(col_weights)
  }
  
  # Spaltenwahl für Baum
  all_pos <- seq_len(aln_len)
  if (!is.null(feature_pos)) {
    feature_pos <- sort(unique(as.integer(feature_pos)))
    feature_pos <- feature_pos[feature_pos >= 1L & feature_pos <= aln_len]
    feature_pos <- intersect(feature_pos, keep_cols_global)
  } else {
    feature_pos <- integer(0L)
  }
  
  if (!length(feature_pos)) {
    # Keine expliziten Features -> Baum auf allen Spalten, die Gap-Filter bestehen
    tree_pos <- keep_cols_global
  } else if (length(feature_pos) >= feature_min_for_subsample) {
    tree_pos <- feature_pos
  } else {
    other_pos <- setdiff(keep_cols_global, feature_pos)
    n_extra <- if (length(other_pos) > 0L) {
      ceiling(extra_col_frac * length(keep_cols_global))
    } else 0L
    if (n_extra > length(other_pos)) n_extra <- length(other_pos)
    extra_pos <- if (n_extra > 0L) sort(sample(other_pos, n_extra, replace = FALSE)) else integer(0L)
    tree_pos  <- sort(unique(c(feature_pos, extra_pos)))
    if (!length(tree_pos)) tree_pos <- keep_cols_global
  }
  
  # Sequenzen für Baum nur über tree_pos
  sub_mat <- aln_mat[, tree_pos, drop = FALSE]
  subseqs <- apply(sub_mat, 1, paste0, collapse = "")
  
  d  <- pwalign::stringDist(Biostrings::DNAStringSet(subseqs), method = "hamming")
  hc <- hclust(d, method = "average")
  
  # Schrittgrenzen wie gehabt, aber auf hc
  if (step_mode == "quantile_merge_heights") {
    merge_sizes <- rev(seq(nrow(hc$merge)))
    heights     <- rev(hc$height)
    merge_df    <- data.frame(height = heights, size = merge_sizes)
    threshold   <- pmax(min_branch_size, min_branch_frac * n)
    filtered    <- merge_df[merge_df$size >= threshold, , drop = FALSE]
    
    probs      <- seq(0, 1, length.out = n_steps + 1)^quantile_skew
    max_height <- max(hc$height)
    
    if (nrow(filtered) < 2L) {
      step_bounds <- c(0, max_height)
    } else {
      raw_bounds <- stats::quantile(filtered$height, probs = probs, names = FALSE)
      min_step_height <- min_step_height_frac * max_height
      keep <- c(TRUE, diff(raw_bounds) > min_step_height)
      step_bounds <- raw_bounds[keep]
      if (length(step_bounds) < 2L) step_bounds <- c(0, max_height)
    }
  } else {
    max_height  <- max(hc$height)
    step_bounds <- seq(0, max_height, length.out = n_steps + 1)
  }
  
  # Vorbereitung für Gewichtung
  col_weights <- numeric(aln_len)
  hit_counts  <- integer(aln_len)
  col_total_contributions <- integer(aln_len)
  
  # Hauptschleife über Schnitte
  for (step in seq_len(length(step_bounds) - 1L)) {
    h_high   <- step_bounds[step + 1L]
    clusters <- cutree(hc, h = h_high)
    groups   <- split(seq_along(clusters), clusters)
    
    edge_hit   <- logical(aln_len)
    edge_score <- numeric(aln_len)
    
    for (idx in groups) {
      m <- length(idx)
      if (m < min_node_size || m > (n - min_node_size)) next
      
      comp_idx <- setdiff(seq_len(n), idx)
      if (!length(comp_idx)) next
      
      g_mat <- base_factor_mat[idx, , drop = FALSE]
      r_mat <- base_factor_mat[comp_idx, , drop = FALSE]
      
      # A C G T nur
      nb <- length(base_levels)
      tab_g <- apply(g_mat, 2, tabulate, nbins = nb)
      tab_r <- apply(r_mat, 2, tabulate, nbins = nb)
      if (!is.matrix(tab_g)) tab_g <- matrix(tab_g, nrow = nb)
      if (!is.matrix(tab_r)) tab_r <- matrix(tab_r, nrow = nb)
      rownames(tab_g) <- base_levels
      rownames(tab_r) <- base_levels
      
      counts_g <- tab_g[c("A","C","G","T"), , drop = FALSE]
      counts_r <- tab_r[c("A","C","G","T"), , drop = FALSE]
      
      depth_g <- colSums(counts_g)
      depth_r <- colSums(counts_r)
      
      ok_depth <- which(depth_g >= min_minor_count & depth_r >= min_minor_count)
      if (!length(ok_depth)) next
      
      # Frequenzen
      freq_g <- sweep(counts_g[, ok_depth, drop = FALSE], 2, pmax(1L, depth_g[ok_depth]), "/")
      freq_r <- sweep(counts_r[, ok_depth, drop = FALSE], 2, pmax(1L, depth_r[ok_depth]), "/")
      
      # Edge-Score pro Spalte: 0.5 * Sum(|p_g - p_r|) in [0,1]
      diff_mat <- abs(freq_g - freq_r)
      score_vec <- 0.5 * colSums(diff_mat)
      
      # Schwelle: starker Kontrast
      # Verbindung zu min_minor_frac: Score >= min_minor_frac
      sel <- ok_depth[score_vec >= min_minor_frac]
      if (!length(sel)) next
      
      # Gap-Filter global
      sel <- intersect(sel, keep_cols_global)
      if (!length(sel)) next
      
      # Score speichern: max über alle Gruppen an diesem Step
      for (j in sel) {
        edge_hit[j]   <- TRUE
        edge_score[j] <- max(edge_score[j], score_vec[match(j, ok_depth)])
      }
    }
    
    sel_cols <- which(edge_hit)
    if (length(sel_cols)) {
      step_weight_base <- 1 / n_steps
      n_informative    <- length(sel_cols)
      for (j in sel_cols) {
        hit_counts[j] <- hit_counts[j] + 1L
        col_total_contributions[j] <- col_total_contributions[j] + n_informative
        # Score skalieren mit Edge-Score
        col_weights[j] <- col_weights[j] + step_weight_base * edge_score[j] / hit_counts[j]
      }
    }
  }
  
  # optionaler globaler Floor wie bisher
  if (apply_global_floor) {
    informative_mask   <- col_weights > 0
    n_informative      <- sum(informative_mask)
    n_uninformative    <- sum(!informative_mask)
    
    if (n_uninformative > 0L && n_informative > 0L) {
      total_informative_weight <- sum(col_weights[informative_mask])
      max_allowed_floor        <- 0.1 * max(col_weights[informative_mask])
      proposed_floor_weight    <- total_informative_weight / n_uninformative
      
      floor_weight <- if (proposed_floor_weight > max_allowed_floor) {
        max_allowed_floor
      } else {
        proposed_floor_weight
      }
      col_weights[!informative_mask] <- floor_weight
    }
  }
  
  names(col_weights) <- seq_len(aln_len)
  attr(col_weights, "hit_counts")        <- hit_counts
  attr(col_weights, "avg_cols_per_step") <-
    ifelse(hit_counts > 0, col_total_contributions / hit_counts, NA_real_)
  attr(col_weights, "step_heights")      <- step_bounds
  attr(col_weights, "tree_positions")    <- tree_pos
  attr(col_weights, "feature_positions") <- feature_pos
  
  col_weights * 1000
}


prepare_one_hot_matrix <- function(fq,
                                   weights,
                                   one_hot_function = one_hot_matrix_conservative) {
  
  Rcpp::cppFunction('
NumericMatrix one_hot_matrix_conservative(CharacterVector alignment,
                                          CharacterVector qualities,
                                          NumericVector weights = NumericVector::create()) {
  int n_seq = alignment.size();
  int aln_len = as<std::string>(alignment[0]).length();

  std::vector<std::string> qual_vec(n_seq);
  for (int i = 0; i < n_seq; ++i) {
    qual_vec[i] = Rcpp::as<std::string>(qualities[i]);
  }

  std::vector<std::map<char, int>> base_counts(aln_len);
  for (int j = 0; j < aln_len; ++j) {
    std::map<char, int>& counts = base_counts[j];
    for (int i = 0; i < n_seq; ++i) {
      char base = alignment[i][j];
      if (base != \'N\') counts[base]++;
    }
  }

  std::vector<int> final_cols;
  std::vector<std::vector<char>> final_bases;
  for (int j = 0; j < aln_len; ++j) {
    std::vector<char> keep;
    for (auto& p : base_counts[j]) keep.push_back(p.first);
    if (!keep.empty()) {
      final_cols.push_back(j);
      final_bases.push_back(keep);
    }
  }

  int num_features = 0;
  for (auto& vec : final_bases) num_features += vec.size();
  NumericMatrix mat(n_seq, num_features);
  CharacterVector colnames(num_features);
  CharacterVector rownames(n_seq);

  int offset = 0;
  for (size_t i = 0; i < final_cols.size(); ++i) {
    int pos = final_cols[i];
    std::map<char, int> base_map;
    for (size_t k = 0; k < final_bases[i].size(); ++k) {
      base_map[final_bases[i][k]] = k;
      colnames[offset + k] = std::string(1, final_bases[i][k]) + "_" + std::to_string(pos);
    }

    for (int r = 0; r < n_seq; ++r) {
      char base = alignment[r][pos];
      double confidence;
      if (base == \'-\') {
        confidence = 0.98;
      } else {
        char qchar = qual_vec[r][pos];
        int qual_score = static_cast<int>(qchar) - 33;
        double prob = std::pow(10.0, -qual_score / 10.0);
        confidence = std::pow(1.0 - prob, 2.0);
      }
      if (base_map.count(base)) {
        double weight = (weights.size() == aln_len) ? weights[pos] : 1.0;
        mat(r, offset + base_map[base]) = confidence * weight;
      }
    }
    offset += final_bases[i].size();
  }

  CharacterVector aln_names;
  bool has_names = false;
  if (alignment.hasAttribute("names")) {
    aln_names = alignment.attr("names");
    has_names = aln_names.size() == n_seq;
  }
  for (int r = 0; r < n_seq; ++r) {
    rownames[r] = has_names ? Rcpp::as<std::string>(aln_names[r]) : std::to_string(r + 1);
  }

  Rcpp::colnames(mat) = colnames;
  Rcpp::rownames(mat) = rownames;
  return mat;
}
')
  
  
  # Extract aligned sequences and quality strings
  aligned_seqs <- as.character(ShortRead::sread(fq))
  qualities <- as.character(quality(fq@quality))
  
  # Generate the one-hot matrix using the provided function
  mat <- one_hot_function(
    alignment = aligned_seqs,
    qualities = qualities,
    weights = weights
  )
  
  # Assign rownames from FASTQ sequence IDs
  rownames(mat) <- as.character(fq@id)
  return(mat)
}


# --- Hybrid distance: combine strengths of dist1 (bridge rescue) and dist2 (IUPAC control) ---
weighted_jaccard_dist <- function(feat, wgts){
  # Inputs: binary one‑hot features (rows=reads, cols=pos:base), column weights
  # Output: R 'dist' that adapts per pair: 
  #  - High overlap → behave like weighted_jaccard_dist2 (fewer bridges → fewer IUPAC)
  #  - Low overlap or thin rows → borrow from weighted_jaccard_dist1 (rescue bridges → less noise)
  
  # --- helpers --------------------------------------------------------------
  mk_dist <- function(M, method_tag = "wj_hybrid"){ n <- nrow(M); lt <- M[lower.tri(M, FALSE)];A
  fin <- is.finite(lt); if(!any(fin)) stop("all distances non-finite"); med <- stats::median(lt[fin]);
  lt[!fin] <- med; lt[lt < 0] <- 0; attr(lt,"Size")<-n; attr(lt,"Diag")<-FALSE;
  attr(lt,"Upper")<-FALSE; attr(lt,"method")<-method_tag; class(lt)<-"dist"; lt }
  v2m <- function(v, n){ M <- matrix(0, n, n); k <- 1L; for(j in 1:(n-1L)) for(i in (j+1L):n){ M[i,j]<-v[k]; M[j,i]<-v[k]; k<-k+1L }; M }
  rank_pow <- function(M, pow){ lt <- M[lower.tri(M, FALSE)]; F <- stats::ecdf(lt[is.finite(lt)]); v <- F(lt)^pow; v2m(v, nrow(M)) }
  
  # --- input prep -----------------------------------------------------------
  feat <- as.matrix(feat); wgts <- as.numeric(wgts); wgts[!is.finite(wgts)|wgts<0] <- 0
  if(!length(wgts) || sum(wgts)==0) wgts <- rep(1, ncol(feat))
  N <- nrow(feat); if(N < 2L) stop("N<2"); J <- ncol(feat)
  s_w <- sqrt(wgts / pmax(1e-12, sum(wgts)))
  
  # --- base components (shared for both worlds) -----------------------------
  if(!exists("dist_from_A_Uc", mode="function")) Rcpp::cppFunction('
    #include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP dist_from_A_Uc(NumericMatrix A, NumericMatrix Uc, double eps=1e-12){
  int n=A.nrow(); R_xlen_t L=(R_xlen_t)n*(n-1)/2, k=0; NumericVector d(L);
  for(int j=0;j<n-1;++j){ for(int i=j+1;i<n;++i){ double den = 1.0 - Uc(i,j); if(!R_finite(den) || den < eps) den = eps;
    double S = A(i,j) / den; if(!R_finite(S)) S = 0.0; d[k++] = 1.0 - S; }}
  d.attr("Size")=n; d.attr("Diag")=false; d.attr("Upper")=false; d.attr("method")="wj_fused"; d.attr("class")="dist"; return d; }')
  
  Xw  <- sweep(feat,     2, s_w, "*")
  Xwc <- sweep(1.0-feat, 2, s_w, "*")
  A   <- tcrossprod(Xw)             # agreement on present features
  Uc  <- tcrossprod(Xwc)            # agreement on complements
  Dfus<- dist_from_A_Uc(A, Uc)      # 1 - A/(1-Uc)
  Mfus<- v2m(as.numeric(Dfus), N); diag(Mfus) <- 0
  
  binZ <- (feat > 0L)
  Pw   <- sweep(binZ*1L, 2, s_w, "*")
  Spos <- tcrossprod(Pw)
  rpw  <- rowSums(Pw^2)             # row load in weighted one‑hot space
  Upw  <- outer(rpw, rpw, `+`) - Spos
  Mpos <- 1 - (Spos / pmax(1e-12, Upw)); diag(Mpos) <- 0
  Mpos[!is.finite(Mpos)] <- stats::median(Mpos[is.finite(Mpos)])
  
  TH   <- cbind(Xw, Xwc)            # signed two‑hot cosine
  TH   <- sweep(TH, 2, colMeans(TH), "-")
  G    <- tcrossprod(TH)
  n2   <- sqrt(pmax(1e-12, rowSums(TH^2)))
  COS  <- G / tcrossprod(n2)
  Msig <- 1 - COS; diag(Msig) <- 0
  Msig[!is.finite(Msig)] <- stats::median(Msig[is.finite(Msig)])
  
  # --- global diagnostics ---------------------------------------------------
  prev    <- pmin(pmax(colMeans(binZ), 1e-12), 1-1e-12)
  entr    <- -(prev*log(prev) + (1-prev)*log1p(-prev)) / log(2)
  lt_fus  <- Mfus[lower.tri(Mfus)]
  ecf     <- stats::ecdf(lt_fus[is.finite(lt_fus)])
  r_med   <- stats::median(ecf(lt_fus)); r_iqr <- stats::IQR(ecf(lt_fus))
  deficit <- r_med * (1 - r_iqr)
  midm    <- mean(prev >= stats::quantile(prev,.25) & prev <= stats::quantile(prev,.75))
  w_pos1  <- pmin(0.5, pmax(0, mean(entr) * deficit * (1 - midm)))
  
  lt_sig  <- Msig[lower.tri(Msig)]
  adv_s   <- mean(pmax(0, lt_sig - lt_fus), na.rm=TRUE)
  adv_f   <- mean(pmax(0, lt_fus - lt_sig), na.rm=TRUE)
  adv_rt  <- if((adv_s+adv_f)>0) adv_s/(adv_s+adv_f) else 0
  nnz     <- colSums(feat>0L)
  s_spr   <- 1 - stats::median(nnz / N)                   # column sparsity proxy
  
  # --- World 1: "bridge rescue" mix (simplified dist1) ---------------------
  w_rem1 <- 1 - w_pos1
  w_sig1 <- w_rem1 * (adv_rt) * 0.5
  w_fus1 <- w_rem1 - w_sig1
  M1     <- w_fus1*Mfus + w_pos1*Mpos + w_sig1*Msig
  if(deficit > 0){ M1 <- rank_pow(M1, 1/(1+deficit)); lt <- M1[lower.tri(M1)]; d0 <- stats::quantile(lt, probs=1/max(2L,length(lt)), na.rm=TRUE); M1[M1<d0] <- d0; diag(M1)<-0 }
  
  # --- World 2: "IUPAC control" mix (dist2 core with gates) ----------------
  tau   <- stats::median(Upw[is.finite(Upw) & Upw>0]); if(!is.finite(tau)||tau<=0) tau <- 1
  Gov   <- Upw / (Upw + tau)                            # overlap gate [0..1]
  Adv   <- pmax(0, Msig - Mfus)
  Advn  <- Adv / pmax(1e-12, abs(Msig) + abs(Mfus))
  lambda<- 0.5 * adv_rt * s_spr
  w_rem2<- 1 - w_pos1
  w_sig2<- w_rem2 * (adv_rt) * 0.5
  w_fus2<- w_rem2 - w_sig2
  Boost <- pmin(w_fus2, lambda * (1 - Gov) * Advn * w_rem2)
  Wsig  <- w_sig2 + Boost
  Wfus  <- w_fus2 - Boost
  Mpair <- Wfus*Mfus + w_pos1*Mpos + Wsig*Msig
  q25   <- stats::quantile(rpw[rpw>0], 0.25, na.rm=TRUE); if(!is.finite(q25)||q25<=0) q25 <- 1
  gi    <- pmin(1, rpw / q25)
  Grow  <- 0.5 + 0.5*sqrt(outer(gi, gi))
  Mbad  <- w_pos1*Mpos + (w_rem2)*Msig
  M2    <- Grow * ( Gov * Mpair + (1 - Gov) * Mbad )
  if(deficit > 0){ M2 <- rank_pow(M2, 1/(1+deficit)); lt <- M2[lower.tri(M2)]; d0 <- stats::quantile(lt, probs=1/max(2L,length(lt)), na.rm=TRUE); M2[M2<d0] <- d0; diag(M2)<-0 }
  
  # --- Hybrid gate: choose world by overlap and thinness --------------------
  Gthin <- 1 - sqrt(outer(gi, gi))                      # thin pairs → 1
  alpha <- pmin(1, pmax(0, 0.5*(1 - Gov) + 0.5*Gthin)) # weight for M1
  Mh    <- alpha*M1 + (1 - alpha)*M2
  
  # Diagnostics out
  options(
    .wj_last_weights = c(
      fused1 = w_fus1, pos1 = w_pos1, signed1 = w_sig1,
      fused2 = mean(Wfus[lower.tri(Wfus)]), pos2 = w_pos1, signed2 = mean(Wsig[lower.tri(Wsig)]),
      alpha_med = stats::median(alpha[lower.tri(alpha)])
    ),
    .wj_last_meta = list(N=N, J=J, sparsity=s_spr, adv_ratio=adv_rt, tau=tau, deficit=deficit)
  )
  
  mk_dist(Mh, "wj_hybrid")
}

# 1.4 distance mix (Jaccard + kmer/UMAP) with per-read coverage
mix_dist_auto_noknn <- function(
    um,
    D_jaccard,
    features,
    floor_kmer    = 0.5,
    lambda_min    = 0.05,
    lambda_max    = 0.25,
    min_feat_read = minpts_sweep.min_depth
){
  U <- as.matrix(um)
  if (!inherits(D_jaccard, "dist")) {
    stop("D_jaccard must be a 'dist' object")
  }
  n <- nrow(U)
  if (is.null(n) || n <= 1L) {
    stop("um must be an N x p matrix with N > 1")
  }
  D_um <- dist(U)
  dj <- as.numeric(D_jaccard)
  du <- as.numeric(D_um)
  med_j <- stats::median(dj[is.finite(dj)])
  med_u <- stats::median(du[is.finite(du)])
  if (!is.finite(med_j) || med_j == 0) med_j <- 1
  if (!is.finite(med_u) || med_u == 0) med_u <- 1
  dj <- dj / med_j
  du <- du / med_u
  F <- as.matrix(features)
  if (nrow(F) != n) {
    stop("features must have same number of rows as um")
  }
  cov_raw  <- rowSums(F > 0L)
  cov_norm <- cov_raw / pmax(1L, ncol(F))
  lam_vec <- lambda_max - (lambda_max - lambda_min) * pmin(1, cov_norm * ncol(F) / pmax(1L, min_feat_read))
  out <- numeric(length(dj))
  k <- 1L
  for (i in 2:n){
    for (j in 1:(i - 1L)){
      lam_j <- max(lambda_min, min(lambda_max, max(lam_vec[i], lam_vec[j])))
      lam_u <- max(floor_kmer, 1 - lam_j)
      s <- lam_j + lam_u
      out[k] <- (lam_j * dj[k] + lam_u * du[k]) / s
      k <- k + 1L
    }
  }
  attr(out, "Size")   <- n
  attr(out, "Diag")   <- FALSE
  attr(out, "Upper")  <- FALSE
  attr(out, "method") <- "mix_jaccard_umap"
  class(out) <- "dist"
  out
}

# 1.4 # -----------------------------------------------------------------------
# Adaptive minPts sweep with dynamic refinement and safeguards
# -----------------------------------------------------------------------------
minpts_sweep2 <- function(
    D,
    grid         = minpts_sweep.grid,
    fq,
    ambig_thr    = ambiguity_by_cluster.ambig_thr,
    min_depth    = ambiguity_by_cluster.min_depth,
    weights      = minpts_sweep.weights,
    extra_minPts = c(7L, 5L, 3L)
){
  distanceMatrixInput    <- D
  minPointsGridValues    <- grid
  fastqObject            <- fq
  ambiguityThreshold     <- ambig_thr
  minimumDepth           <- min_depth
  objectiveWeightVector  <- weights
  extraMinPointsSequence <- as.integer(extra_minPts)
  
  absoluteAmbiguityGuard   <- 10L
  fractionalAmbiguityGuard <- 0.10
  hardClusterReadCap       <- 100L
  
  # Normalize and sanitize initial minPts grid
  if (is.function(minPointsGridValues)) {
    minPointsGridValues <- minPointsGridValues()
  }
  minPointsGridValues <- as.integer(minPointsGridValues)
  minPointsGridValues <- minPointsGridValues[is.finite(minPointsGridValues) & minPointsGridValues > 0L]
  minPointsGridValues <- sort(unique(minPointsGridValues))
  
  # Basic sequence information
  dnaStringSetAll        <- Biostrings::DNAStringSet(ShortRead::sread(fastqObject))
  totalReadCount         <- length(dnaStringSetAll)
  bigClusterSizeLimit    <- 0.10 * totalReadCount
  maxAllowedMinPts       <- totalReadCount -1
  
  # Rule (1): drop grid points larger than 2x number of reads
  minPointsGridValues <- minPointsGridValues[minPointsGridValues <= maxAllowedMinPts]
  
  # If nothing survives, fall back to extra_minPts ("small extension points")
  if (!length(minPointsGridValues)) {
    extraMinPointsSequence <- sort(unique(extraMinPointsSequence[extraMinPointsSequence >= 3L &
                                                                   extraMinPointsSequence <= maxAllowedMinPts]))
    minPointsGridValues <- extraMinPointsSequence
  }
  
  # If still empty, return an empty scored data.frame
  if (!length(minPointsGridValues)) {
    return(data.frame(
      minPts = integer(0), n_clusters = integer(0), max_size = integer(0),
      noise_frac = numeric(0), prob_mean = numeric(0), stability_mean = numeric(0),
      ambig_bases = numeric(0), ambig_bases_raw = numeric(0),
      s_noise = numeric(0), s_prob = numeric(0), s_stab = numeric(0), s_ambig = numeric(0),
      s_parsimony = numeric(0), objective = numeric(0)
    ))
  }
  
  initialMinPointsGrid <- minPointsGridValues
  
  # Global column filtering
  consensusMatrixAll      <- consensus_matrix_ACGTN_gap(dnaStringSetAll)
  gapFractionAllColumns   <- as.numeric(consensusMatrixAll["-", ] / pmax(1L, colSums(consensusMatrixAll)))
  keepColumnIndicesGlobal <- which(gapFractionAllColumns <= 0.90)
  minimumDepthAdjusted    <- max(10L, as.integer(minimumDepth))
  
  # Count IUPAC-like ambiguous positions in a DNAStringSet subset
  countIupacLocal <- function(dnaSubset){
    count_iupac_cols(
      dna_set   = dnaSubset,
      keep_idx  = keepColumnIndicesGlobal,
      depth_min = minimumDepthAdjusted,
      thr       = ambiguityThreshold,
      min_count = 3L
    )
  }
  
  # Evaluate clustering quality metrics for a single minPts choice
  evaluateMinPointsValue <- function(minPointsValue){
    hdbscanFitObject          <- dbscan::hdbscan(distanceMatrixInput, minPts = minPointsValue, cluster_selection_epsilon = 0.05)
    clusterLabelVector        <- as.integer(hdbscanFitObject$cluster)
    clusterSizeVector         <- as.integer(table(clusterLabelVector[clusterLabelVector > 0L]))
    clusterCountValue         <- length(clusterSizeVector)
    maximumClusterSize        <- if (length(clusterSizeVector)) max(clusterSizeVector) else 0L
    noiseFraction             <- mean(clusterLabelVector == 0L)
    meanMembershipProbability <- mean(hdbscanFitObject$membership_prob,                        na.rm = TRUE)
    meanStabilityProbability  <- mean(hdbscanFitObject$membership_prob[hdbscanFitObject$cluster > 0L], na.rm = TRUE)
    inClusterReadCount        <- sum(clusterLabelVector > 0L)
    
    if (clusterCountValue > 0L){
      rawAmbiguousBaseCount      <- 0L
      weightedAmbiguousBaseCount <- 0
      clusterIdentifierVector    <- sort(setdiff(unique(clusterLabelVector), 0L))
      skipReadBudget             <- as.integer(round(fractionalAmbiguityGuard * maximumClusterSize))
      skippedReadCount           <- 0L
      for (clusterIdentifier in clusterIdentifierVector){
        readIndexVector  <- which(clusterLabelVector == clusterIdentifier)
        clusterReadCount <- length(readIndexVector)
        if (clusterReadCount < absoluteAmbiguityGuard) {
          next
        }
        if (clusterReadCount < hardClusterReadCap && skippedReadCount + clusterReadCount <= skipReadBudget) {
          skippedReadCount <- skippedReadCount + clusterReadCount
          next
        }
        ambiguousCountCluster      <- countIupacLocal(dnaStringSetAll[readIndexVector])
        rawAmbiguousBaseCount      <- rawAmbiguousBaseCount + ambiguousCountCluster
        weightedAmbiguousBaseCount <- weightedAmbiguousBaseCount + (clusterReadCount / totalReadCount) * ambiguousCountCluster
      }
    } else {
      rawAmbiguousBaseCount      <- countIupacLocal(dnaStringSetAll)
      weightedAmbiguousBaseCount <- rawAmbiguousBaseCount
    }
    
    data.frame(
      minPts            = as.integer(minPointsValue),
      n_clusters        = clusterCountValue,
      max_size          = maximumClusterSize,
      noise_frac        = noiseFraction,
      prob_mean         = meanMembershipProbability,
      stability_mean    = meanStabilityProbability,
      ambig_bases       = weightedAmbiguousBaseCount,
      ambig_bases_raw   = rawAmbiguousBaseCount,
      in_cluster_reads  = inClusterReadCount,
      stringsAsFactors  = FALSE,
      check.names       = FALSE
    )
  }
  
  # Safe row binding with harmonized column sets
  bindRowsSafe <- function(mainDataFrame, extraDataFrame){
    extraDataFrame <- as.data.frame(extraDataFrame, stringsAsFactors = FALSE, check.names = FALSE)
    
    extraOnlyNames <- setdiff(names(extraDataFrame), names(mainDataFrame))
    for (columnName in extraOnlyNames) {
      mainDataFrame[[columnName]] <- NA
    }
    
    mainOnlyNames <- setdiff(names(mainDataFrame), names(extraDataFrame))
    for (columnName in mainOnlyNames) {
      extraDataFrame[[columnName]] <- NA
    }
    
    extraDataFrame <- extraDataFrame[, names(mainDataFrame), drop = FALSE]
    rbind(mainDataFrame, extraDataFrame)
  }
  
  # Linear rescaling into [0,1] with NA-safe handling
  scaleToUnitInterval <- function(numericVectorInput){
    finiteMask <- is.finite(numericVectorInput)
    if (!any(finiteMask)) return(rep(1, length(numericVectorInput)))
    rangeValues <- range(numericVectorInput[finiteMask])
    if (diff(rangeValues) == 0) {
      scaledVector <- rep(1, length(numericVectorInput))
    } else {
      scaledVector <- (numericVectorInput - rangeValues[1]) / diff(rangeValues)
    }
    scaledVector[!finiteMask] <- 1
    scaledVector
  }
  
  # Score all grid rows with noise, stability, ambiguity and parsimony terms
  scoreResultDataFrame <- function(resultDataFrame){
    weightVectorDefault <- c(noise = 0, prob = 0, stability = 0, ambig = 0, parsimony = 0)
    weightVectorDefault[names(objectiveWeightVector)] <- as.numeric(objectiveWeightVector)
    weightVectorDenominator <- sum(abs(weightVectorDefault))
    if (!is.finite(weightVectorDenominator) || weightVectorDenominator == 0) weightVectorDenominator <- 1
    weightVectorNormalized <- weightVectorDefault / weightVectorDenominator
    
    resultDataFrame$s_noise <- scaleToUnitInterval(1 - resultDataFrame$noise_frac)
    resultDataFrame$s_prob  <- scaleToUnitInterval(replace(resultDataFrame$prob_mean,      !is.finite(resultDataFrame$prob_mean),      0))
    resultDataFrame$s_stab  <- scaleToUnitInterval(replace(resultDataFrame$stability_mean, !is.finite(resultDataFrame$stability_mean), 0))
    
    ambiguousWeighted <- replace(resultDataFrame$ambig_bases,     !is.finite(resultDataFrame$ambig_bases),     0)
    ambiguousRaw      <- replace(resultDataFrame$ambig_bases_raw, !is.finite(resultDataFrame$ambig_bases_raw), 0)
    ambiguousWeighted <- pmax(0, as.numeric(ambiguousWeighted))
    ambiguousRaw      <- pmax(0, as.numeric(ambiguousRaw))
    
    scaleAmbiguity <- function(ambiguityVector){
      if (!length(ambiguityVector)) return(rep(0, length(ambiguityVector)))
      minValue <- min(ambiguityVector, na.rm = TRUE)
      maxValue <- max(ambiguityVector, na.rm = TRUE)
      if (!is.finite(maxValue) || maxValue <= minValue) return(rep(0, length(ambiguityVector)))
      scaledVector <- (ambiguityVector - minValue) / (maxValue - minValue)
      scaledVector[scaledVector < 0] <- 0
      scaledVector[scaledVector > 1] <- 1
      scaledVector
    }
    
    ambiguityScaledWeighted <- scaleAmbiguity(ambiguousWeighted)
    ambiguityScaledRaw      <- scaleAmbiguity(ambiguousRaw)
    
    combinedAmbiguity       <- 0.5 * ambiguityScaledWeighted + 0.5 * ambiguityScaledRaw
    ambiguityScore          <- 1 - (combinedAmbiguity ^ 0.6)
    moderateMask            <- combinedAmbiguity >= 0.50 & combinedAmbiguity <= 0.60
    ambiguityScore[moderateMask] <- pmin(ambiguityScore[moderateMask], 0.80)
    mildMask                <- combinedAmbiguity >= 0.35 & combinedAmbiguity < 0.50
    ambiguityScore[mildMask] <- pmin(ambiguityScore[mildMask], 0.90)
    resultDataFrame$s_ambig <- ambiguityScore
    
    # Use precomputed in-cluster read counts (no second hdbscan call)
    inClusterReadCounts <- resultDataFrame$in_cluster_reads
    
    maxShareVector  <- ifelse(
      inClusterReadCounts > 0,
      resultDataFrame$max_size / pmax(inClusterReadCounts, resultDataFrame$minPts * 10L),
      0
    )
    minPointsFactor <- 1 / log1p(resultDataFrame$minPts)
    parsimonyRaw    <- maxShareVector * minPointsFactor
    resultDataFrame$s_parsimony <- scaleToUnitInterval(parsimonyRaw)
    
    resultDataFrame$objective <- rowSums(cbind(
      resultDataFrame$s_noise     * weightVectorNormalized["noise"],
      resultDataFrame$s_prob      * weightVectorNormalized["prob"],
      resultDataFrame$s_stab      * weightVectorNormalized["stability"],
      resultDataFrame$s_ambig     * weightVectorNormalized["ambig"],
      resultDataFrame$s_parsimony * weightVectorNormalized["parsimony"]
    ), na.rm = TRUE)
    
    resultDataFrame$objective[resultDataFrame$n_clusters == 0L] <-
      0.5 * resultDataFrame$objective[resultDataFrame$n_clusters == 0L]
    
    resultDataFrame
  }
  
  # Evaluate grid adaptively with early stopping
  resultDataFrame   <- NULL
  evaluatedCount    <- 0L
  zeroClusterCount  <- 0L
  
  L          <- length(minPointsGridValues)
  halfGrid   <- ceiling(L / 2)
  allIndices <- seq_len(L)
  
  # Determine up to three middle indices
  if (L >= 3L) {
    centerIndex <- floor((L + 1L) / 2L)
    middleIndices <- sort(unique(pmax(1L, pmin(L, centerIndex + c(-1L, 0L, 1L)))))
  } else {
    middleIndices <- allIndices
  }
  
  remainingIndices <- setdiff(allIndices, middleIndices)
  
  appendRow <- function(rowDf){
    if (is.null(resultDataFrame)) {
      resultDataFrame <<- rowDf
    } else {
      resultDataFrame <<- bindRowsSafe(resultDataFrame, rowDf)
    }
    evaluatedCount <<- evaluatedCount + 1L
    if (rowDf$n_clusters[1] == 0L) {
      zeroClusterCount <<- zeroClusterCount + 1L
    }
  }
  
  # Evaluate middle grid points first
  for (idx in middleIndices) {
    mp  <- minPointsGridValues[idx]
    row <- evaluateMinPointsValue(mp)
    appendRow(row)
  }
  
  # Rule (2): if three middle points all look perfect, stop
  if (length(middleIndices) >= 3L) {
    midMinPts     <- minPointsGridValues[middleIndices]
    midMask       <- resultDataFrame$minPts %in% midMinPts
    midRows       <- resultDataFrame[midMask, , drop = FALSE]
    if (nrow(midRows) >= 3L) {
      sameClusterCount <- length(unique(midRows$n_clusters)) == 1L
      positiveClusters <- midRows$n_clusters[1] > 0L
      allIupacZero     <- all(midRows$ambig_bases_raw == 0L)
      if (sameClusterCount && positiveClusters && allIupacZero) {
        resultDataFrame <- scoreResultDataFrame(resultDataFrame)
        return(resultDataFrame[order(-resultDataFrame$objective, resultDataFrame$minPts),
                               c("minPts", "n_clusters", "max_size", "noise_frac", "prob_mean", "stability_mean",
                                 "ambig_bases", "ambig_bases_raw", "s_noise", "s_prob", "s_stab", "s_ambig",
                                 "s_parsimony", "objective")])
      }
    }
  }
  
  # Evaluate remaining grid points with early stop on pure noise
  for (idx in remainingIndices) {
    mp  <- minPointsGridValues[idx]
    row <- evaluateMinPointsValue(mp)
    appendRow(row)
    
    # Rule (3): if at least half of the grid has been checked and all seen so far are pure noise, stop
    if (evaluatedCount >= halfGrid && zeroClusterCount == evaluatedCount) {
      break
    }
  }
  
  # First scoring on the evaluated grid
  resultDataFrame <- scoreResultDataFrame(resultDataFrame)
  
  # Check whether the initial grid already supports a sufficiently large cluster
  hasInitialBigCluster <- any(
    resultDataFrame$minPts %in% initialMinPointsGrid &
      resultDataFrame$n_clusters > 0L &
      resultDataFrame$max_size >= bigClusterSizeLimit
  )
  
  # Safeguard: if the best configuration already has clusters and zero raw IUPAC count,
  # then dynamic refinement to smaller minPts is skipped.
  bestInitialObjectiveIndex    <- which.max(resultDataFrame$objective)
  bestInitialHasClusters       <- resultDataFrame$n_clusters[bestInitialObjectiveIndex] > 0L
  bestInitialRawIupacCount     <- resultDataFrame$ambig_bases_raw[bestInitialObjectiveIndex]
  skipRefinementDueToCleanBest <- bestInitialHasClusters &&
    is.finite(bestInitialRawIupacCount) && bestInitialRawIupacCount == 0
  
  positiveClusterDataFrame <- resultDataFrame[resultDataFrame$n_clusters > 0L, , drop = FALSE]
  
  if (nrow(positiveClusterDataFrame) > 0L && !skipRefinementDueToCleanBest) {
    # Case 1: at least one configuration with non-zero clusters and refinement is allowed
    bestPositiveIndex         <- which.max(positiveClusterDataFrame$objective)
    bestPositiveMinPoints     <- positiveClusterDataFrame$minPts[bestPositiveIndex]
    bestPositiveAmbiguity     <- positiveClusterDataFrame$ambig_bases[bestPositiveIndex]
    smallestPositiveMinPoints <- min(positiveClusterDataFrame$minPts)
    isBestAtLeftEdge          <- (bestPositiveMinPoints == smallestPositiveMinPoints)
    
    ambiguityValues        <- positiveClusterDataFrame$ambig_bases
    nonZeroAmbiguityValues <- ambiguityValues[ambiguityValues > 0 & is.finite(ambiguityValues)]
    if (length(nonZeroAmbiguityValues)) {
      medianAmbiguity <- stats::median(nonZeroAmbiguityValues)
      isAmbiguityHigh <- is.finite(bestPositiveAmbiguity) && bestPositiveAmbiguity >= medianAmbiguity
    } else {
      isAmbiguityHigh <- FALSE
    }
    
    clusterCountValues <- positiveClusterDataFrame$n_clusters
    medianClusterCount <- stats::median(clusterCountValues)
    isClusterCountLow  <- positiveClusterDataFrame$n_clusters[bestPositiveIndex] <= medianClusterCount
    
    if (isBestAtLeftEdge && isAmbiguityHigh && isClusterCountLow) {
      # Refine grid downward using configurable extra_minPts
      extraMinPointsSequence <- sort(unique(extraMinPointsSequence[extraMinPointsSequence > 0L]))
      
      if (length(extraMinPointsSequence) > 1L) {
        firstStageExtraMinPoints  <- extraMinPointsSequence[-length(extraMinPointsSequence)]
        secondStageExtraMinPoints <- extraMinPointsSequence[length(extraMinPointsSequence)]
      } else {
        firstStageExtraMinPoints  <- extraMinPointsSequence
        secondStageExtraMinPoints <- integer(0L)
      }
      
      firstStageExtraMinPoints <- firstStageExtraMinPoints[
        firstStageExtraMinPoints < bestPositiveMinPoints & firstStageExtraMinPoints >= 3L
      ]
      firstStageExtraMinPoints <- setdiff(firstStageExtraMinPoints, resultDataFrame$minPts)
      
      if (length(firstStageExtraMinPoints)) {
        firstStageRowsList       <- lapply(firstStageExtraMinPoints, evaluateMinPointsValue)
        firstStageExtraDataFrame <- do.call(rbind, firstStageRowsList)
        resultDataFrame          <- bindRowsSafe(resultDataFrame, firstStageExtraDataFrame)
        resultDataFrame          <- scoreResultDataFrame(resultDataFrame)
        
        positiveClusterDataFrame2 <- resultDataFrame[resultDataFrame$n_clusters > 0L, , drop = FALSE]
        
        if (nrow(positiveClusterDataFrame2) > 0L && length(secondStageExtraMinPoints) == 1L) {
          bestStageTwoIndex     <- which.max(positiveClusterDataFrame2$objective)
          bestStageTwoAmbiguity <- positiveClusterDataFrame2$ambig_bases[bestStageTwoIndex]
          
          if (is.finite(bestPositiveAmbiguity) && is.finite(bestStageTwoAmbiguity) &&
              bestStageTwoAmbiguity >= bestPositiveAmbiguity) {
            
            secondStageCandidateValue <- as.integer(secondStageExtraMinPoints)
            if (!secondStageCandidateValue %in% resultDataFrame$minPts &&
                secondStageCandidateValue >= 3L &&
                secondStageCandidateValue < bestPositiveMinPoints) {
              secondStageExtraDataFrame <- evaluateMinPointsValue(secondStageCandidateValue)
              resultDataFrame           <- bindRowsSafe(resultDataFrame, secondStageExtraDataFrame)
              resultDataFrame           <- scoreResultDataFrame(resultDataFrame)
            }
          }
        }
      }
    }
    
  } else if (nrow(positiveClusterDataFrame) == 0L) {
    # Case 2: all evaluated grid values currently lead to pure noise
    smallestCurrentMinPoints <- min(resultDataFrame$minPts)
    if (smallestCurrentMinPoints > 3L) {
      lowerNewMinPoints      <- max(3L, floor(smallestCurrentMinPoints / 2L))
      newMinPointsGridValues <- seq(lowerNewMinPoints, smallestCurrentMinPoints - 1L)
      newMinPointsGridValues <- setdiff(newMinPointsGridValues, resultDataFrame$minPts)
      if (length(newMinPointsGridValues)) {
        newRowsList       <- lapply(newMinPointsGridValues, evaluateMinPointsValue)
        newExtraDataFrame <- do.call(rbind, newRowsList)
        resultDataFrame   <- bindRowsSafe(resultDataFrame, newExtraDataFrame)
        resultDataFrame   <- scoreResultDataFrame(resultDataFrame)
      }
    }
  }
  
  # Safeguard: newly added smaller minPts values must not destroy large clusters
  hasInitialClusterAboveTenPercent <- any(
    resultDataFrame$minPts %in% initialMinPointsGrid &
      resultDataFrame$n_clusters > 0L &
      resultDataFrame$max_size >= bigClusterSizeLimit
  )
  hasInitialClusterAboveSix <- any(
    resultDataFrame$minPts %in% initialMinPointsGrid &
      resultDataFrame$n_clusters > 0L &
      resultDataFrame$max_size >= 6L
  )
  
  if (hasInitialClusterAboveTenPercent || hasInitialClusterAboveSix) {
    isNewSmallerMinPoints <- resultDataFrame$minPts < min(initialMinPointsGrid) &
      !(resultDataFrame$minPts %in% initialMinPointsGrid)
    
    isSmallClusterRelative <- hasInitialClusterAboveTenPercent &
      isNewSmallerMinPoints & resultDataFrame$max_size < bigClusterSizeLimit
    
    isSmallClusterAbsolute <- hasInitialClusterAboveSix &
      isNewSmallerMinPoints & resultDataFrame$max_size < 6L
    
    isSmallClusterAtNewMinPts <- isSmallClusterRelative | isSmallClusterAbsolute
    
    if (any(isSmallClusterAtNewMinPts)) {
      minimumObjectiveValue <- min(resultDataFrame$objective, na.rm = TRUE)
      resultDataFrame$objective[isSmallClusterAtNewMinPts] <- minimumObjectiveValue - 1
    }
  }
  
  resultDataFrame[order(-resultDataFrame$objective, resultDataFrame$minPts),
                  c("minPts", "n_clusters", "max_size", "noise_frac", "prob_mean", "stability_mean",
                    "ambig_bases", "ambig_bases_raw", "s_noise", "s_prob", "s_stab", "s_ambig",
                    "s_parsimony", "objective")]
}


# 1.2 IUPAC counter for a DNA subset with global column filter
count_iupac_cols <- function(dna_set,
                             keep_idx,
                             depth_min = minpts_sweep.min_depth,
                             thr = ambiguity_by_cluster.ambig_thr,
                             min_count = 3L){
  if (length(dna_set) == 0L) return(0L)
  cm  <- consensus_matrix_ACGTN_gap(dna_set)
  M4  <- cm[c("A","C","G","T"), , drop = FALSE]
  dep <- colSums(M4)
  kk  <- intersect(keep_idx, which(dep >= depth_min))
  if (!length(kk)) return(0L)
  Mkk <- M4[, kk, drop = FALSE]
  dkk <- colSums(Mkk)
  Fkk <- sweep(Mkk, 2, pmax(1L, dkk), "/")
  sum(vapply(seq_along(kk), function(i){
    vv   <- sort(Fkk[, i],   decreasing = TRUE)
    cnts <- sort(Mkk[, i],   decreasing = TRUE)
    if (length(vv) < 2L) return(0L)
    if (cnts[2L] < min_count) return(0L)
    as.integer(vv[2L] >= thr)
  }, integer(1L)))
}

plot_umap_clusters_3d_from_results <- function(results,
                                               um,
                                               file_prefix = "umap3d_clusters",
                                               out_dir = ".",
                                               point_size = 4) {
  stopifnot(is.list(results), "final" %in% names(results))
  
  df <- results$final
  stopifnot(all(c("sequence", "cluster_id") %in% colnames(df)))
  
  ## 1) extract coordinates ----------------------------------------------------
  if (is.matrix(um) || is.data.frame(um)) {
    umap_coords <- as.matrix(um)
  } else if (!is.null(um$layout)) {
    umap_coords <- as.matrix(um$layout)
  } else if (!is.null(um$embedding)) {
    umap_coords <- as.matrix(um$embedding)
  } else {
    stop("Cannot extract coordinates from 'um'.")
  }
  
  if (ncol(umap_coords) < 3)
    stop("UMAP object must have at least 3 dimensions for a 3D plot.")
  
  if (is.null(rownames(umap_coords)))
    stop("UMAP coordinates need rownames matching results$final$sequence.")
  
  ## 2) align IDs -------------------------------------------------------------
  cl_vec <- setNames(df$cluster_id, df$sequence)
  common_ids <- intersect(rownames(umap_coords), names(cl_vec))
  if (!length(common_ids))
    stop("No overlap between rownames(um) and results$final$sequence.")
  
  common_ids <- sort(common_ids)
  coords3d <- umap_coords[common_ids, 1:3, drop = FALSE]
  clusters <- cl_vec[common_ids]
  
  plot_df <- data.frame(
    UMAP1      = coords3d[, 1],
    UMAP2      = coords3d[, 2],
    UMAP3      = coords3d[, 3],
    cluster_id = as.integer(clusters),
    sequence   = common_ids,
    stringsAsFactors = FALSE
  )
  
  ## 3) color mapping, fully opaque -------------------------------------------
  suppressPackageStartupMessages({
    library(plotly)
    library(htmlwidgets)
    library(RColorBrewer)
  })
  
  cl_ids <- sort(unique(plot_df$cluster_id))
  n_cl   <- length(cl_ids)
  
  base_pal <- brewer.pal(min(max(n_cl - ("0" %in% cl_ids)), 12), "Set3")
  
  col_map <- setNames(rep(NA_character_, length(cl_ids)), cl_ids)
  if (0 %in% cl_ids) col_map[as.character(0)] <- "#B0B0B0"  # noise = solid grey
  nonzero <- setdiff(cl_ids, 0)
  if (length(nonzero) > 0) {
    col_map[as.character(nonzero)] <- rep(base_pal, length.out = length(nonzero))
  }
  
  plot_df$cluster_factor <- factor(plot_df$cluster_id,
                                   levels = cl_ids,
                                   ordered = TRUE)
  
  colors_vec <- as.vector(col_map[levels(plot_df$cluster_factor)])
  
  ## 4) 3D plot with larger fonts / legend ------------------------------------
  p3d <- plot_ly(
    data  = plot_df,
    x     = ~UMAP1,
    y     = ~UMAP2,
    z     = ~UMAP3,
    type  = "scatter3d",
    mode  = "markers",
    color = ~cluster_factor,
    colors = colors_vec,
    text  = ~sequence,
    hoverinfo = "text",
    marker = list(
      size = point_size,                 # also controls legend marker size
      line = list(color = "black", width = 0.7),
      opacity = 1                        # fully opaque, no alpha
    )
  ) %>%
    layout(
      title = list(
        text = "UMAP 3D projection with cluster assignments",
        font = list(size = 26)
      ),
      scene = list(
        xaxis = list(
          title    = list(text = "UMAP 1", font = list(size = 20)),
          tickfont = list(size = 16)
        ),
        yaxis = list(
          title    = list(text = "UMAP 2", font = list(size = 20)),
          tickfont = list(size = 16)
        ),
        zaxis = list(
          title    = list(text = "UMAP 3", font = list(size = 20)),
          tickfont = list(size = 16)
        )
      ),
      legend = list(
        title = list(text = "Cluster ID", font = list(size = 20)),
        font  = list(size = 18),          # legend text ~4× typical default
        itemsizing = "constant"
      ),
      font = list(size = 16)             # global default for other text
    )
  
  ## 5) save html --------------------------------------------------------------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  html_file <- file.path(out_dir, paste0(file_prefix, "_UMAP3D_clusters.html"))
  saveWidget(as_widget(p3d), html_file, selfcontained = TRUE)
  message("Saved 3D UMAP cluster plot to: ", html_file)
  
  invisible(list(plot = p3d, file = html_file))
}

