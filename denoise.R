if (!"this.path" %in% installed.packages()) install.packages("this.path")
library(this.path)

options(install.packages.compile.from.source = "never")
Sys.setenv(R_COMPILE_AND_INSTALL_PACKAGES = "never")

source(file = paste0(this.path::here(),"/scripts/startup.R"))
source(file = paste0(this.path::here(),"/scripts/build_tools.R"))

check_and_fix_rcpp_build(verbose = TRUE, prefer = "ucrt64") #set path manually: rtools_dir = "C:/Rtools45"

source(file = paste0(this.path::here(),"/scripts/align_fastq.R"))
source(file = paste0(this.path::here(),"/scripts/filter_nanopore_batch.R"))
source(file = paste0(this.path::here(),"/scripts/png_gallery_to_pdf_paginated.R"))
source(file = paste0(this.path::here(),"/scripts/homopolymer.R"))
source(file = paste0(this.path::here(),"/scripts/clustering.R"))
source(file = paste0(this.path::here(),"/scripts/post_merge_v3.R"))
source(file = paste0(this.path::here(),"/scripts/consensus_first_pass_v4.R"))

required <- c(
  "ShortRead", "Biostrings", "Rcpp", "ggplot2", "doParallel", "foreach", "dbscan", "reshape2",
  "magrittr", "Matrix", "dplyr", "stringr", "plotly", "htmlwidgets",
  "RColorBrewer", "this.path", "uwot", "tidyr", "png", "pkgbuild")

#in case MKL is installed, default to single thread
Sys.setenv(
  MKL_THREADING_LAYER = "SEQUENTIAL",
  MKL_NUM_THREADS     = "1",
  OMP_NUM_THREADS     = "1",
  MKL_DYNAMIC         = "FALSE",
  KMP_AFFINITY        = "none",
  OPENBLAS_NUM_THREADS= "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

ensure_required_packages(required)

path_actual <- this.path::here()
setwd(path_actual)

source(file = paste0(this.path::here(),"/config.txt"))

if (mode=="fastq")
{
fastq_files <- list.files(path = sample_path, pattern = paste0(input_suffix, ".fastq"), full.names = TRUE)
fastq_files <- validate_nanopore_filter_inputs(fastq_files, n_thread)[[1]]

batch_summary <- filter_nanopore_batch(fastq_files, plot = TRUE,
                                       n_cores = n_thread,
                                       filter_nanopore_batch.suffix_lowq,
                                       filter_nanopore_batch.max_removal_fraction,
                                       filter_nanopore_batch.mad_threshold)

png_gallery_to_pdf_paginated(folder = sample_path, filter_nanopore_batch.output_quality_pdf)
}

if (mode=="fasta")
{
  
  fasta_to_fastq_phred <- function(fasta_files, phred_score = 20) {
    qchar <- intToUtf8(phred_score + 33)
    fasta_files %>%
      lapply(function(f) {
        fa <- readFasta(f)
        seqs <- sread(fa) %>% gsub("-", "", .) %>% DNAStringSet
        quals <- BStringSet(
          vapply(width(seqs), function(n) paste(rep(qchar, n), collapse = ""), character(1))
        )
        fq <- ShortReadQ(sread = seqs, quality = quals, id = ShortRead::id(fa))
        out <- sub("\\.fa(sta)?$", paste0(filter_nanopore_batch.suffix_lowq,".fastq"), f, ignore.case = TRUE)
        if (file.exists(out)) file.remove(out)
        writeFastq(fq, out, compress = FALSE)
      })
  }
  
  fasta_files <- list.files(path = sample_path, pattern = paste0(input_suffix, ".fasta"), full.names = TRUE)
  invisible(fasta_to_fastq_phred(fasta_files, phred_score = dummy_qscores_phred))
  
  # fasta_files <- list.files(path = sample_path, pattern = "ONT.*refined.fasta", full.names = TRUE)
  # fasta_to_fastq_phred(fasta_files, phred_score = 22)
  # 
  # fasta_files <- list.files(path = sample_path, pattern = "PB.*refined.fasta", full.names = TRUE)
  # fasta_to_fastq_phred(fasta_files, phred_score = 40)
  
}


##########################FASTQ ALIGNMENT
fastq_files <- list.files(path = sample_path, pattern = paste0(filter_nanopore_batch.suffix_lowq,".fastq$") , full.names = TRUE)

closeAllConnections()
cl <- makeCluster(n_thread) #make sure to have 1.5 GB per 500 reads (per thread)
registerDoParallel(cl)

results_align <- foreach(i = seq_along(fastq_files), .inorder = FALSE, .packages = c("ShortRead", "Biostrings", "magrittr")) %dopar% {
                     align_fastq_with_mafft(
                       fastq_file = fastq_files[i],
                       path_mafft = align_fastq_with_mafft.my_mafft,
                       output_suffix = align_fastq_with_mafft.output_suffix_align,
                       param_threads = 1, #fast alignment method --retree scales poorly
                       min_avg_qual = align_fastq_with_mafft.min_avg_qual,
                       gap_quality_char = align_fastq_with_mafft.gap_quality_char,
                       overwrite = align_fastq_with_mafft.overwrite_align,
                       verbose = TRUE)
                   }
stopCluster(cl); closeAllConnections()


##########################HOMOPOLYMER REMOVAL
fastq_files <- list.files(path = sample_path, pattern = paste0(align_fastq_with_mafft.output_suffix_align,".fastq$"), full.names = TRUE)

cl <- makeCluster(n_thread) #mafft inside foreach also uses up multiple threads
registerDoParallel(cl)

results_nohomo <- foreach(i = seq_along(fastq_files), .inorder = FALSE, .packages = c("ShortRead", "Rcpp")) %dopar% {
  filter_homopolymer_alignment(
    fastq_path = fastq_files[i], fastq_dir = sample_path,
    min_len = filter_homopolymer_alignment.min_len,
    homopolymer_threshold = filter_homopolymer_alignment.homopolymer_threshold,
    gap_fraction_cutoff = filter_homopolymer_alignment.gap_fraction_cutoff,
    output_suffix = filter_homopolymer_alignment.output_suffix,
    di_rep = 6,
    tri_rep = 6,
    use_di_tri = FALSE)
}

stopCluster(cl); closeAllConnections()


##########################NOISE MASKING
fastq_files <- list.files(path = sample_path, pattern = paste0(filter_homopolymer_alignment.output_suffix,".fastq$"), full.names = TRUE)

cl <- parallel::makeCluster(n_thread) #mafft inside foreach also uses up multiple threads
registerDoParallel(cl)


results_cluster <- foreach(i = seq_along(fastq_files),
                           .packages = c("ShortRead", "Biostrings", "dbscan", "Rcpp", "reshape2", "magrittr", "plotly", "htmlwidgets", "RColorBrewer", "uwot"),
                           .inorder = TRUE) %dopar% 
{                             
  file_path <- fastq_files[i]
  fq <- readFastq(file_path)
  
  tryCatch({
    if (!file.exists(file_path)) {
      return(paste("File does not exist: ", file_path))
    }
    
    seqs <- as.character(sread(fq))
    if (length(seqs) < 10) {
      return(paste("Too few sequences (<10) in: ", file_path))
    }

    set.seed(1)
    # extract features
    extr_features <- detect_rare_minor_dual(
      fq            = fq,
      window        = detect_rare_minor_dual.window,
      max_gap_frac  = detect_rare_minor_dual.max_gap_frac,
      min_frac      = detect_rare_minor_dual.min_frac,
      min_reads     = detect_rare_minor_dual.min_reads,
      fdr           = detect_rare_minor_dual.fdr,
      use_median_bg = detect_rare_minor_dual.use_median_bg
    )
    
    feat_pos <- extr_features$features %>% colnames %>% gsub("pos(.*):.*","\\1",.) %>% as.integer

    w_col <- hierarchical_column_weights(
      aligned_seqs = Biostrings::DNAStringSet(ShortRead::sread(fq)),
      feature_pos  = feat_pos,
      min_node_size = ambiguity_by_cluster.min_depth,
      quantile_skew = 5,
      include_gaps    = FALSE,
      max_gap_fraction = detect_rare_minor_dual.max_gap_frac,
      n_steps = 50,
      min_step_height_frac = 0.01,
      min_branch_frac = 0.1,
      min_branch_size = round(width(fq)[1]/100)
    )

    one_hot_matr <- prepare_one_hot_matrix(fq, w_col)
    set.seed(1)
    um <- uwot::umap2(one_hot_matr, n_neighbors = min(20, round(length(fq)/2)), n_components = 5, min_dist = 0)

    if (is.null(extr_features$features) || ncol(extr_features$features) == 0L) {
      ## if no features = dummy matrix
      n        <- length(fq)
      feat_mat <- matrix(0L, nrow = n, ncol = 1L)
      D_feat   <- dist(matrix(0, nrow = n, ncol = 1L))
    } else {
      ## weighted Jaccard distance from features
      feat_mat <- extr_features$features
      w <- extr_features$weights
      if (!length(w) || length(w) != ncol(feat_mat) || !any(is.finite(w))) {
        w <- rep(1, ncol(feat_mat))
      }
      D_feat <- weighted_jaccard_dist(feat_mat, w)
    }
    
    # merge distance metrics
    dmix <- mix_dist_auto_noknn(
      um              = um,
      D_jaccard       = D_feat,
      features        = feat_mat,
      floor_kmer      = 0.2,
      lambda_min      = 0.5,
      lambda_max      = 0.8,
      min_feat_read   = ambiguity_by_cluster.min_depth
    )

    # minPts Sweep
    my_minpts <- minpts_sweep2(
      D          = dmix,
      grid       = minpts_sweep.grid,
      fq         = fq,
      ambig_thr  = ambiguity_by_cluster.ambig_thr,
      min_depth  = ambiguity_by_cluster.min_depth,
      weights    = minpts_sweep.weights,
      extra_minPts = c(7L, 5L, 4L, 3L)
    )

    sel_minPts <- as.integer(my_minpts$minPts[1L])
    
    #  HDBSCAN
    set.seed(1)
    fit  <- dbscan::hdbscan(dmix, minPts = sel_minPts, cluster_selection_epsilon = 0.05)
    labs <- as.integer(fit$cluster)
    
    # Outlier-Scores 
    seq_ids <- as.character(ShortRead::id(fq))
    os <- fit$outlier_scores
    if (is.null(os) || !length(os)) {
      os <- rep(0, length(labs))
    }
    
    is_outlier <- as.integer(labs == 0L | os >= detect_rare_minor_dual.max_gap_frac)
    
    # if all reads are noise = outlier, set all reads to cluster 1 and remove outlier flag
    if (all(labs == 0L) && all(is_outlier == 1L)) {
      labs       <- rep(1L, length(labs))
      is_outlier <- rep(0L, length(labs))
    }
    
    # Ambiguity per cluster:
    dna_all     <- Biostrings::DNAStringSet(ShortRead::sread(fq))
    cm_all      <- consensus_matrix_ACGTN_gap(dna_all)
    keep_global <- which(as.numeric(cm_all["-", ] / pmax(1L, colSums(cm_all))) <= detect_rare_minor_dual.max_gap_frac)
    cl_ids      <- sort(setdiff(unique(labs), 0L))
    ambig_list  <- vector("list", length(cl_ids) + 1L)
    
    ambig_list[[1L]] <- data.frame(
      cluster     = 0L,
      n_reads     = sum(labs == 0L),
      ambig_sites = NA_integer_,
      stringsAsFactors = FALSE
    )
    
    if (length(cl_ids)) {
      for (i in seq_along(cl_ids)) {
        cl    <- cl_ids[i]
        idx_c <- which(labs == cl)
        amb   <- count_iupac_cols(
          dna_all[idx_c],
          keep_idx  = keep_global,
          depth_min = ambiguity_by_cluster.min_depth,
          thr       = ambiguity_by_cluster.ambig_thr,
          min_count = detect_rare_minor_dual.min_reads
        )
        ambig_list[[i + 1L]] <- data.frame(
          cluster     = cl,
          n_reads     = length(idx_c),
          ambig_sites = amb,
          stringsAsFactors = FALSE
        )
      }
    }
    ambig_base <- do.call(rbind, ambig_list)
    
    # final output
    final_df <- data.frame(
      sequence   = seq_ids,
      cluster_id = labs,
      is_outlier = is_outlier,
      stringsAsFactors = FALSE
    )
    
    results <- list(
      minpts     = my_minpts,
      base       = labs,
      ambig_base = ambig_base,
      final      = final_df
    )
    
    if(plotting_html_3d)
    {
      um <- uwot::umap2(one_hot_matr, n_neighbors = min(20, round(length(fq)/2)), n_components = 3, min_dist = 0)
      plot_umap_clusters_3d_from_results(
        results     = results,
        um          = um,
        file_prefix = file_path %>% gsub(paste0(".*/|_",input_suffix,".*"),"",.),
        out_dir     = paste0(sample_path, "plots")
      )
    }
    
    return(results$final)
      
  }, error = function(e) {
    return(paste("Error processing file: ", file_path, " | ", conditionMessage(e)))
  })
}


names(results_cluster) <- tools::file_path_sans_ext(basename(fastq_files)) %>% gsub(filter_homopolymer_alignment.output_suffix,"",.)

stopCluster(cl); closeAllConnections() #ignore warnings about closed connections

#lengths(results_cluster) %>% table

saveRDS(results_cluster, paste0(sample_path, "results_cluster.rds"))

to_delete <- list.files(
  path = sample_path,
  pattern = "_CL",
  full.names = TRUE
)

if (length(to_delete) > 0) {
  suppressWarnings(invisible(file.remove(to_delete)))
}

to_delete <- list.files(
  path = sample_path,
  pattern = "_final",
  full.names = TRUE
)

if (length(to_delete) > 0) {
  suppressWarnings(invisible(file.remove(to_delete)))
}

cl <- makeCluster(n_thread)
registerDoParallel(cl)

# Parallel consensus generation
results_consensus <- foreach(i = seq_along(results_cluster), 
                             .packages = c("Biostrings", "ShortRead","magrittr", "dplyr", "stringr"), 
                             .inorder = TRUE) %dopar% 
  {
    cluster_name <- names(results_cluster)[i]
    cluster_data <- results_cluster[[i]]
    colnames(cluster_data) <- c("sequence","cluster_id","is_outlier")
    
    # Basic checks
    if (is.null(cluster_data) || !is.data.frame(cluster_data) || nrow(cluster_data) < 6) {
      # Fix skipped case
      return(paste("Skipping", cluster_name, ": too few sequences or no data."))
    }
    
    tryCatch({
      generate_consensus_with_export(
        results_cluster_i = setNames(list(cluster_data), cluster_name),
        fastq_dir =  sample_path,
        min_cluster_size = generate_consensus_with_export.min_cluster_size,
        output_suffix = generate_consensus_with_export.output_suffix,
        iupac_min_count = generate_consensus_with_export.iupac_min_count,
        iupac_min_prop = generate_consensus_with_export.iupac_min_prop
      )
      return(cluster_name)  # Return processed name as success
    }, error = function(e) {
      return(paste("Error in", cluster_name, ":", conditionMessage(e)))
    })
  }

stopCluster(cl)
closeAllConnections()

# Clean results: names of all successfully processed clusters
results_consensus <- Filter(Negate(is.null), results_consensus)


#post-merging
fasta_files <- list.files(sample_path, pattern = paste0(generate_consensus_with_export.output_suffix,"\\.fasta$"), full.names = TRUE)
fasta_files <- fasta_files[( basename(fasta_files) %>% sub(paste0("_", input_suffix,".*"),"",.) ) %in% ( results_consensus %>% unlist %>% sub(paste0("_", input_suffix,".*"),"",.) )]
fasta_files <- fasta_files[ !grepl("_final_merged\\.fasta$", basename(fasta_files)) ]

cl <- makeCluster(n_thread)
registerDoParallel(cl)

results_merged <- foreach(i = seq_along(fasta_files),
                          .packages = c("Biostrings", "stringr", "tools", "tibble", "Rcpp"),
                          .inorder = FALSE) %dopar% 
  {
    
    input_file <- fasta_files[i]
    
    # Defensive try block
    tryCatch({
      if (!file.exists(input_file) || file.info(input_file)$size == 0) {
        return(paste("Skipping", basename(input_file), ": file does not exist or is empty."))
      }
      
      # Optionally also check number of sequences inside
      seqs <- Biostrings::readDNAStringSet(input_file)
      if (length(seqs) == 0) {
        return(paste("Skipping", basename(input_file), ": no sequences in FASTA."))
      }

      # Run merging logic
      res <- merge_and_regenerate_consensus(
        input_fasta_file = input_file,
        iupac_min_prop = merge_and_regenerate_consensus.iupac_min_prop,
        min_cluster_size = merge_and_regenerate_consensus.min_cluster_size,
        iupac_min_count = merge_and_regenerate_consensus.iupac_min_count,
        max_pairwise_distance = merge_and_regenerate_consensus.max_pairwise_distance
      )
      
      return(res)
    }, error = function(e) {
      return(paste("Error in", basename(input_file), ":", conditionMessage(e)))
    })
  }

stopCluster(cl)
closeAllConnections()

# Extract usable results only (list outputs, not character error messages)
results_merged_clean <- Filter(Negate(is.character), results_merged)

# Optional: inspect messages
error_messages <- Filter(is.character, results_merged)
cat(paste(unlist(error_messages), collapse = "\n"))

summary_df <- collect_final_merged_consensus(input_dir = sample_path, input_suffix = input_suffix)
print(summary_df)

temp_files <- list.files(path=sample_path, "_final_merged.fasta", full.names = TRUE)
all_cluster <- NULL

for (i in seq_along(temp_files))
{
  temp <- temp_files[i] %>% readDNAStringSet
  names(temp) <- names(temp) %>% paste0(temp_files[i] %>% basename %>% gsub(paste0("_", input_suffix,".*"),"",.), "__", .)
  all_cluster <- append(all_cluster, temp)
}

temp_matches <- match(summary_df[,c(1,2)] %>% apply(.,1,paste0,collapse="__"), names(all_cluster) %>% gsub("_IMD.*|_IAD.*","",.))

final_out <- cbind(summary_df, seq=all_cluster[temp_matches])

DNAStringSet(final_out[,"seq"]  %>% gsub("-","",.)) %>% setNames(rownames(final_out)) %>% writeXStringSet(., file=paste0(sample_path, "results.fasta"))


