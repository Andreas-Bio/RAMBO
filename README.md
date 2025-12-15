**This is the repository of: RAMBO: Resolving Amplicons in Mixed Samples for Accurate DNA Barcoding with Oxford Nanopore (doi: to_be_done)**


# RAMBO
Resolving Amplicons in Mixed Samples for Accurate DNA Barcoding with Oxford Nanopore

## Introduction and Intended Use

The RAMBO workflow, started by calling the denoise.R script, is an R-based workflow designed for DNA barcode sequencing analysis (e.g. mitochondrial COI, rRNA genes). It focuses on datasets with high-coverage Nanopore amplicon reads, where the goal is to recover high-fidelity consensus sequences. The pipeline works best if each sample contains reads from one dominant single specimen, but can deal with other cases too.

Key requirements and recommendations include:

1) Minimum reads per sample: Ensure at least ~20 sequences per sample for meaningful extraction of the dominant contig. Additional contigs, each, also require a minimum read count (~10) to be represented in the output files as separate clusters. Fewer reads may result in the pipeline skipping clustering for that sample or producing mixed clusters of secondary contigs.

2) Maximum gap levels: Alignment columns with more than 90% gaps are removed. If an input file contains highly divergent taxa (i.e. multiple families, each equally represented), the alignment can become excessively gapped and most columns may be discarded, leaving no usable alignment for downstream clustering.

3) Memory and threads: Multi-threading is supported. Allocate approximately 2–4 GB of RAM per thread for typical run sizes; as a rule of thumb, about 1.5 GB per 500 Nanopore reads per thread is recommended. Using too many threads without sufficient memory can lead to slowdowns or crashes, as memory usage does not scale linearly. Settings should therefore be adapted to the available hardware. It is strongly recommended to downsample samples with 10,000 reads or more.

4) Data preparation: Input reads should be primer-trimmed, demultiplexed, and length-filtered for plausibility before running the pipeline. Extremely low-quality reads or non-target-length sequences should be removed upstream for best results.


## Setup and Initialization

Supported Platforms: The pipeline is written in R and has been tested on Windows, but is also competible with macOS and Linux systems. It requires a working R installation (R >4.5.0). 

Required R Packages: The pipeline relies on both CRAN and Bioconductor packages. Key dependencies include:
•	Bioconductor packages: ShortRead, Biostrings, and pwalign for reading FASTQ/FASTA data and sequence alignment operations.
•	CRAN packages: Rcpp (for high-performance C++ routines), parallel/doParallel and foreach (for parallel processing loops), dbscan (HDBSCAN clustering algorithm), and various tidyverse/utility libraries (ggplot2, reshape2, dplyr, stringr, tibble, etc.) for data manipulation and plotting. Visualization libraries like plotly, htmlwidgets, and color palettes from RColorBrewer are included for optional outputs. Matrix and vector packages (Matrix, magrittr) are also used. All required package names are defined in the script’s required list.
•	BiocManager: Ensure BiocManager is installed to fetch Bioconductor packages. The pipeline provides a helper that will install any missing packages automatically. On first run, the script will attempt to install ShortRead, Biostrings, pwalign, and all other required packages if they are not already available. This may take some time.

Note that R always requires / in paths regardless of the operating system, while the command prompt/console requires either / or \\\, depending on your OS.

**MAFFT** (multiple sequence alignment program) is invoked by the pipeline for aligning reads. Download MAFFT for your OS and note its executable path. For example, on Windows the path might be C:/Program Files/mafft/mafft.bat, whereas on Linux/macOS you can install via package managers and use the output of which mafft (e.g. /usr/local/bin/mafft). In the configuration (see below), you'll specify the full path to the MAFFT executable. Test by running path/to/mafft --help in a terminal.

**Make**: Required for compiling packages that include native code. This is typically available on Unix systems, if not install the package build-essential. Windows users must install [RTools](https://cran.r-project.org/bin/windows/Rtools/). The pipeline includes a check_and_fix_rcpp_build() function that tries to configure your system’s C++ toolchain automatically. You may see messages about checking for make or gcc – if it reports they are missing, install the appropriate tools and re-run.

Obtain the pipeline scripts: Download or clone the repository containing denoise.R and its companion R script files and folders. Place them in a working directory. All processing is done locally, no internet access is required beyond package installation.

---

### Running the Pipeline

Before you start:
Create a directory and copy all input sequencing files (FASTA or FASTQ) into it. All files must share a common filename suffix, which serves as an anchor for identifying valid input files. By default, this suffix is `_refined`, so filenames should follow the pattern: `sample_name_12345_refined.fastq` Avoid spaces and special characters in file names, use only alphanumeric characters and underscores. The chosen suffix must be specified in the configuration file and is used by the pipeline to unambiguously detect input files within the directory. Make sure to have full read and write access in this folder.

Although the main script installs all R packages automatically, it is recommended to pre-install the required packages, as this step might require user intervention if an installations fails.

Open a new R session and run:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)}

BiocManager::install(
  c("ShortRead", "Biostrings", "pwalign", "IRanges"),
  ask = FALSE,
  update = FALSE,
  quiet = TRUE)

install.packages(
  c("Rcpp", "ggplot2", "doParallel", "foreach", "dbscan", "reshape2",
    "magrittr", "Matrix", "dplyr", "stringr", "plotly", "htmlwidgets",
    "RColorBrewer", "this.path", "uwot", "tidyr", "png", "stringr",
"pkgbuild"),  quiet = TRUE)

```

All packages have been pre-installed and a **renv** folder for the Windows OS can be downloaded from the associated Zenodo depository (doi: to_be_done).

---

### Set up configuration

Open the `config.txt` file and define all parameters. Most parameters never require adjustment, while some should be routinely checked for every run (marked with **!**).
The R code in this file will be sourced by the main script.

**General Parameters**

- **!** `mode`  
  Input file mode. Must be either `"fastq"` or `"fasta"`.  
  Use `"fastq"` if quality scores are available (recommended for Oxford Nanopore data).  
  Use `"fasta"` for any data without quality scores. In this case, dummy quality values are assigned internally.

- **!** `input_suffix`  
  Filename suffix (excluding file extension) shared by all input files.  
  For example, if files are named `sample1_refined.fastq` and `sample2_refined.fastq`, set: `input_suffix = "refined"`
  This suffix is used to identify and group valid input files.

- **!** `sample_path`  
  Directory containing all input files.  
  Make this an absolute path, if possible. A trailing slash is required. `sample_path = "/data/barcode_run1/"`

- `n_thread`  
  Number of CPU threads to use.  
  Default is one less than the total number of available cores and capped at 16.  
  Increasing this speeds up processing of multiple samples in parallel, but requires approximately 2–4 GB RAM per thread for ~500 reads per sample.
  
**FASTQ Quality Filtering Parameters**
Only applicable if `mode = "fastq"`. In `fasta` mode, synthetic quality scores are used.

- `filter_nanopore_batch.max_removal_fraction`  
  Maximum fraction of reads removed based on quality.  
  Default `0.1` allows removal of the worst 10 percent of reads.  
  Lower values are more conservative; higher values allow stronger filtering.

- `filter_nanopore_batch.mad_threshold`  
  Quality cutoff expressed as median absolute deviation (MAD).  
  Default `3` provides a balance between sensitivity and stringency for Nanopore R10 data.

- `filter_nanopore_batch.output_quality_pdf`  
  Output filename for the quality summary PDF.  
  Default is `quality_summary.pdf` in the sample directory.

Other parameters in this block (e.g. `suffix_lowq`, `dummy_qscores_phred`) are internal and should not be modified.


**Alignment Parameters**

- `align_fastq_with_mafft.my_mafft`  
  Full path to the MAFFT executable.  
  Examples:
  ```
  "C:/Users/Name/mafft-win/mafft.bat"
  "/usr/bin/mafft"
  ```

- **!** `align_fastq_with_mafft.min_avg_qual`  
  Minimum mean Phred quality required to retain a read for alignment.  
  Default is `16`, tuned for ONT R10.4 FAST basecalling.  
  Higher values increase stringency; lower values retain more reads at the cost of higher error rates.
  Adjust this value depending on the quality of your data (e.g. up to 30 for PacBio). Please remember that the mean quality per read should be calculated by transforming the log into absolute error values and transforing the mean back to log Phred values.

- `align_fastq_with_mafft.overwrite_align`  
  Logical flag controlling whether existing alignment files are overwritten.  
  Default `TRUE` forces realignment.  
  Set to `FALSE` to resume interrupted runs, but remove any incomplete (e.g. 0-byte) files manually.

Advanced MAFFT parameters (e.g. `param_retree`, `param_threads`) are set internally.  
The pipeline parallelizes at the sample level and uses one MAFFT thread per sample.

**Homopolymer Filtering Parameters**

- `filter_homopolymer_alignment.min_len`  
  Minimum homopolymer length triggering column removal.  
  Default `5`.

- `filter_homopolymer_alignment.homopolymer_threshold`  
  Fraction of sequences required to show a homopolymer for column removal.  
  Default `0.33`.

- `filter_homopolymer_alignment.gap_fraction_cutoff`  
  Gap fraction threshold for removing adjacent alignment columns.  
  Default `0.5`.

Output files from this stage use the suffix `_nohomo.fastq`.

**Clustering and Variant Detection Parameters**

These parameters control detection of rare variants, ambiguity assessment, and clustering behavior.  
Defaults are tuned for barcode-scale amplicon data and generally should not be modified unless the data strongly deviate from these assumptions.

Local Background and Sliding Statistics

- `sliding_stat.window`  
  Half-window size (in alignment columns) used to smooth local background rates of non-major bases.  
  Larger values provide stronger smoothing and reduce sensitivity to local fluctuations; smaller values increase local sensitivity but may amplify noise.

- `sliding_stat.fun`  
  Aggregation function used for background smoothing.  
  Default `median` provides robustness against outliers and local sequencing artifacts.

---

Optional Diagnostic Plotting

- `plotting_html_3d`  
  Logical flag controlling generation of interactive 3D UMAP plots.  
  Default `FALSE`. Note that UMAP represents only part of the composite distance metric; clusters in UMAP space may not exactly reflect final clustering results. This parameter increases processing time considerably.

---

Ambiguity Evaluation Parameters

- `ambiguity_by_cluster.ambig_thr`  
  Frequency threshold for calling ambiguity at a column.  
  A column is considered ambiguous if two or more bases each reach at least this fraction of reads.

- `ambiguity_by_cluster.min_depth`  
  Minimum total A, C, G, T depth required at a column for ambiguity evaluation.  
  Columns below this depth are ignored to prevent unstable ambiguity calls.

---

Rare Minor Variant Detection

These parameters control detection of low-frequency variants above local background noise.

- `detect_rare_minor_dual.max_gap_frac`  
  Columns with a gap fraction above this value are excluded from variant detection.  
  Prevents highly gapped alignment regions from influencing rare variant inference.

- `detect_rare_minor_dual.min_frac`  
  Minimum fractional excess above local background required for a minor base to be considered.

- `detect_rare_minor_dual.min_reads`  
  Absolute minimum number of reads supporting a minor base.  
  Prevents single-read artifacts from being treated as true variants.

- `detect_rare_minor_dual.fdr`  
  Benjamini–Hochberg false discovery rate threshold applied across candidate minor sites.

- `detect_rare_minor_dual.use_median_bg`  
  Logical flag controlling whether the median (robust) or mean background rate is used.  
  Default `TRUE` is recommended for noisy long-read data.

---

HDBSCAN Parameter Optimization

- `minpts_sweep.grid`  
  Candidate values for the HDBSCAN `minPts` parameter.  
  The pipeline evaluates this grid and selects the optimal value based on a multi-criterion objective.

- `minpts_sweep.weights`  
  Weights applied to different optimization criteria, including noise fraction, cluster probability, stability, ambiguity, and parsimony.  
  Higher weight on ambiguity penalizes clusters with excessive within-cluster heterogeneity.

---

Initial Consensus Generation Parameters

These parameters control consensus generation for initial clusters.

- `generate_consensus_with_export.min_cluster_size`  
  Minimum number of sequences required to generate a consensus.  
  Clusters smaller than this threshold are ignored to avoid low-confidence consensus sequences.

- `generate_consensus_with_export.output_suffix`  
  Suffix appended to FASTA files containing consensus sequences from the initial clustering stage.

- `generate_consensus_with_export.iupac_min_count`  
  Minimum number of reads supporting a minor base for inclusion as an IUPAC ambiguity code.

- `generate_consensus_with_export.iupac_min_prop`  
  Minimum fraction of reads supporting a minor base for IUPAC encoding.  
  Used in combination with `iupac_min_count` to suppress spurious ambiguities.

---

Cluster Refinement and Consensus Merging Parameters

These parameters control merging of highly similar clusters and regeneration of consensus sequences.

- `merge_and_regenerate_consensus.min_cluster_size`  
  Minimum number of sequences required to consider a cluster or merged group for consensus regeneration.

- `merge_and_regenerate_consensus.iupac_min_count`  
  Minimum number of supporting reads required for IUPAC ambiguity codes during merged consensus generation.

- `merge_and_regenerate_consensus.iupac_min_prop`  
  Minimum fraction of sequences supporting a minor base for IUPAC inclusion after merging.

- `merge_and_regenerate_consensus.max_pairwise_distance`  
  Maximum allowed pairwise distance between clusters to permit merging.  
  A value of `0` disables merging unless consensus sequences are identical.

---

To start the pipeline, open a new R session to ensure a clean environment. The main driver script can then be executed either interactively from within R or non-interactively from the system shell.

From within R or RStudio:
```
source("path/to/denoise.R")
```

From a terminal or command prompt:
```
Rscript path/to/denoise.R
```

Both methods run the same pipeline. The choice depends on whether you prefer interactive execution or command-line batch processing.

It is recommended to run the script from a directory where all input and output folders are properly structured, and to avoid workspace contamination from previous runs.
The folder specified in the config file should be empty, except for the designated input files in fastq or fasta format.

**Test case**: Unzip the file test.zip into the same directory as denoise.R and run the pipeline with no changes.

---

## Main Processing Steps

The denoise.R pipeline consists of a series of stages that sequentially transform raw sequencing reads into final high-confidence consensus sequences. Below is an overview of the main processing steps, in the order they are executed:

1.	**Initial Quality Filtering** (Optional for FASTQ mode): If `mode = "fastq"`, the pipeline first filters out low-quality reads from each input FASTQ file. It calculates the average Phred quality of each read and removes reads below the `min_avg_qual` threshold. It also uses a MAD-based filter to remove up to a set fraction of the worst reads. By default, at most 10% of reads are dropped based on quality metrics. Reads failing these criteria are excluded from downstream analysis to improve clustering accuracy. A summary PDF (quality_summary.pdf) is generated to show quality distributions before and after filtering. (In `mode = "fasta"`, this step is skipped and all reads are retained and given dummy quality scores.) The output of this step is a cleaned FASTQ for each sample, named `<sample>_lowQremoved.fastq`.
  
2.	**Read Alignment** (MAFFT): All remaining reads for each sample are aligned using MAFFT to produce multiple sequence alignments. This is done sample-by-sample in parallel. Before alignment, extremely short files or empty files are skipped (ensuring we don’t try to align if no reads passed filtering). The MAFFT alignment is run with one thread per alignment (the pipeline automatically handles threading at the sample level). Sequence identifiers are preserved and quality scores are mapped onto the alignment (in gaps, a placeholder quality like "!" is inserted). If any reads had non-unique IDs or other issues, the script will error out (ensuring unique IDs by trimming at first whitespace). Low-quality reads may also be filtered here: by default, if after the initial filter a sample still has a few extremely low-quality reads, those are dropped if their average Q is below the set threshold. Reads counts are logged, and if after filtering a sample has ≤10 reads, the alignment step for that sample is aborted with a warning (too few reads to align reliably). The alignment output for each sample is a FASTQ file named `<sample>_aligned.fastq`, containing all aligned reads.

3.	**Homopolymer Removal**: Long homopolymers in Nanopore data can create artifactual alignment columns (extra inserted gaps or erroneous stretches). In this step, the aligned FASTQ is scanned for columns dominated by homopolymer runs. Any alignment column where a significant fraction of reads (≥33% by default) have a run of ≥5 identical bases is removed. Additionally, columns immediately adjacent to a removed homopolymer column are checked: if they consist of >50% gaps, they are also removed as likely spurious. This process effectively “clips out” the regions in the alignment that are not reliable due to systematic sequencing errors, while preserving genuine sequence differences. The result is a refined alignment FASTQ for each sample without those error-prone columns. These outputs are named `<sample>_aligned_nohomo.fastq` (using the default _nohomo suffix). Each sequence in this file is shorter (some columns removed) but still aligned to one another. This stage greatly reduces noise for downstream variant detection.

4.	**Feature Extraction and Clustering**: This is the core step where true sequence variants are distinguished from sequencing noise. For each sample’s filtered alignment:
- The pipeline identifies rare variant features among the reads. It scans each alignment column for minor alleles (non-major bases) that appear above a background error rate and exceed a minimum count. This uses a sliding window to estimate local error rates and a statistical test (binomial with FDR correction) to find positions where a minor base is overrepresented compared to noise. The result is a set of variant features (specific base changes at specific positions) likely representing real polymorphisms.
- It then encodes these features into a distance matrix. Two complementary approaches are combined:
A.	Dimensionality reduction (UMAP): Each read’s sequence is one-hot encoded and projected into a low-dimensional space (5 dimensions by default) using UMAP.
B.	Feature Jaccard distance: A weighted Jaccard distance is computed between reads based on the presence/absence of the detected rare variant features. This emphasizes differences at informative positions.
These distances are mixed into a single composite distance matrix. Essentially, reads are compared by overall sequence similarity (UMAP space) and by specific rare variant patterns, to improve clustering accuracy.
- Next, the pipeline performs a minPts (min cluster size) sweep for HDBSCAN. It tries a range of candidate `minPts` values (e.g. 9, 12, 15, ..., 60) and evaluates a clustering quality score for each. It selects an optimal `minPts` that balances retaining true clusters vs. labeling noise, using an objective weighting of cluster stability, noise, and parsimony. This adaptive approach means you don’t need to manually choose a cluster size parameter. The algorithm finds one suited to your data.
- HDBSCAN clustering: Using the chosen minPts, HDBSCAN (a density-based clustering algorithm) is run on the distance matrix to assign reads to clusters. HDBSCAN will label some reads as outliers (cluster 0) if they do not confidently belong to any cluster. The pipeline records cluster labels for each read. It also retrieves HDBSCAN’s outlier scores and marks any reads with high outlier score (above a threshold, by default 0.90) as outliers as well. So, cluster “0” represents noise or singletons that were not clustered. A special rule is applied: if all reads ended up as outliers (cluster 0), the pipeline will force them into a single cluster rather than leave everything as noise (this prevents total failure on small samples: if nothing clustered, assume one cluster).
- After clustering, the script evaluates ambiguity within each cluster. It goes through each cluster and counts how many sites would be considered polymorphic (ambiguous) within that cluster, based on the `ambig_thr` (e.g. >=25% frequency for a second base). This yields a count of ambiguous sites per cluster which is later reported.
- The outcome of this step for each sample is a set of initial clusters with their member reads. The pipeline saves an internal R data object (results_cluster.rds) containing a data frame for each sample: listing each read, its cluster assignment, and whether it was marked as outlier. For each sample, you can think of this as the initial grouping of reads by putative true sequence. Typically, reads from the same true variant will cluster together, and obvious sequencing errors either join the correct cluster (if minor) or become outliers if they are unique.

5.	**Consensus Sequence Generation** (per sample): For each sample’s clusters, the pipeline now creates consensus sequences to represent each unique variant. This step uses the multiple sequence alignment with homopolymers included:
- Ignoring cluster 0 (noise/outliers) and any clusters smaller than the minimum size (by default, clusters with <6 inlier reads are skipped to avoid low-confidence consensus). It focuses on the valid clusters.
- For each cluster, all reads labeled as inliers (`is_outlier = 0`) are multiple-aligned (they already are from previous steps), and a consensus is calculated. The consensus calling considers the quality scores and the IUPAC thresholds you set: if a minor allele meets the count and frequency criteria, an IUPAC ambiguity code is placed at that position in the consensus sequence. Otherwise, the majority base is used. Any reads in the cluster that were marked as outliers (`is_outlier = 1`) are excluded from consensus calculation to avoid polluting the consensus with potential errors. However, the pipeline keeps track of them (appending “_OUTLIER” to their IDs in an output file) for reference.
- The pipeline also computes various statistics for the cluster:
A. Intra-cluster Average Divergence (IAD): the mean percentage difference of each read to the consensus.
B. Intra-cluster Median Divergence (IMD): the median pairwise difference between reads in the cluster (or equivalently, median distance between reads, which is a measure of cluster tightness).
C. Number of IUPAC codes (NIUPAC): how many ambiguous bases (R, Y, N, etc.) the consensus contains, which is indicative of within-cluster heterogeneity.
D. Number of sequences (NSEQ): the count of reads in the cluster (inliers).
E. Consensus length (LEN): length of the consensus sequence (excluding gaps).
-	The outputs of this step include:
A. Per-cluster FASTA files: For each cluster, a FASTA file with all its member sequences is written (with original read IDs, outliers marked) for auditing. These are named like `<sample>_CL<clusterID>.fasta`. For example, `SampleA_CL1.fasta` contains all reads of cluster 1 from SampleA. (These files are mostly for troubleshooting or if you want to inspect the raw reads of a cluster; they can be large and are optional.)
B. Consensus FASTA file (per sample): A FASTA file summarizing the consensus sequences of all clusters in that sample, named `<sample>_final.fasta` (using the default _final suffix). This is the primary output of the first pass. Each consensus sequence header contains an identifier with the cluster number and the metrics described above. For example, a header might look like: `>CL1_IAD0.0025_IMD0.00_NIUPAC1_NSEQ50_LEN658`
This denotes cluster 1 has IAD = 0.0025 (0.25% avg divergence), IMD = 0.00% (median divergence ~0), 1 IUPAC ambiguity in the consensus, 50 reads, and length 658 bp. These IAD/IMD tags allow quick assessment of cluster quality: low values mean the cluster is very tight (good), whereas higher values might indicate the cluster has substantial internal variation (or possibly mis-clustering). Please note that the IAD and IMD values are typically around your expected error rate and should ideally not be larger than your expected error rate times 2.

After this step, for each sample you have a set of consensus sequences representing the distinct barcode sequences found in that sample’s reads (except that some might still be very similar to each other, which the next step will address).

6.	**Consensus Merging** : The pipeline includes a post-processing step to merge very similar consensus sequences within each sample, which helps if the algorithm accidentally over-split a true sequence into two clusters (for instance, due to a borderline difference or a few errors). This step takes all consensus sequences from a sample’s `_final.fasta` and compares them pairwise:
- It will align consensus sequences and calculate the divergence between them (considering ambiguity codes as matching appropriately). If two consensus sequences differ by less than or equal to the `max_pairwise_distance threshold` (default allows at most 1 mismatch or ambiguous site difference), they are candidates for merging.
- The merging algorithm checks if new ambiguous positions would be introduced by merging clusters. Essentially, this ensures we don’t merge clusters that would create a highly ambiguous consensus (implying they were distinct).
- This process can merge multiple clusters iteratively. For example, if Cluster 1 and Cluster 3 differ only by one base, they merge into a new cluster “CL1_3”. The consensus stats (IAD, IMD, etc.) are recalculated for the merged cluster. The naming convention in the output reflects merged clusters by concatenating their IDs: e.g., `CL1_3_IMD0.50_MDC0.25_NIUPAC2_NSEQ100_LEN658` could denote a merged cluster from original CL1 and CL3, now with updated metrics. Here you’ll notice a new MDC tag in the header, Median Distance to Consensus, which is the median divergence of reads to the new merged consensus (as opposed to IMD which was median between reads). Both IMD and MDC are provided for merged clusters.
- The outputs from this stage are:
A. Merged consensus FASTA (per sample): named `<sample>_final_merged.fasta`, containing the final consensus sequences after any intra-sample merging. If no merges occurred for that sample, this file will be essentially identical to the `_final.fasta` or might contain the same sequences with slightly adjusted names/metrics. If merges did occur, you will see fewer consensus sequences (merged ones) with combined cluster IDs in their names.
B. Merged cluster read files: similar to earlier per-cluster files, the pipeline outputs `<sample>_CL<i>_<j>_merged.fasta` files containing all reads from merged cluster i_j. These are mainly for completeness and debugging, to see which original clusters went into a merged group.

8.	**Final Output Compilation**: In the last step, the pipeline collates all final consensus sequences from all samples into a single multi-FASTA. It reads each sample’s `_final_merged.fasta` and combines them. Each sequence name in this combined FASTA is prefixed with the sample name to keep them unique. For example, if sample “SampleA” had a cluster CL1_IAD..., in the combined file it might appear as SampleA__CL1_IAD.... The combined FASTA is saved as `results.fasta` in the sample_path directory. This results.fasta is convenient for downstream analyses (e.g. BLASTing the consensuses, creating phylogenetic trees, etc.), as it contains one entry per final cluster across all samples.

---

To recap the file outputs you will see after a successful run (for each sample, assuming default suffixes):

•	`<sample>_lowQremoved.fastq`: Cleaned reads after initial quality filter (FASTQ mode only).

•	`<sample>_aligned.fastq`: MAFFT-aligned reads with quality info.

•	`<sample>_aligned_nohomo.fastq`: Aligned reads after homopolymer/gap artifact removal.

•	`<sample>_CL*.fasta`: (Multiple files) Raw reads belonging to each cluster (produced in consensus step).

•	`<sample>_final.fasta`: Consensus sequences for each cluster before merging (with IAD, IMD, etc. in headers).

•	`<sample>_final_merged.fasta`: Final consensus sequences after merging (with updated IMD, MDC in headers).

•	`results.fasta`: Combined all-sample final consensus sequences (each header prefixed by sample and including metrics).

•	Quality and report files: quality_summary.pdf (if FASTQ filtering done) and possibly a log or summary CSV of cluster metrics.

These are the primary outputs. The intermediate FASTQ files (_lowQremoved, _aligned, _aligned_nohomo) can be large; you may delete them if not needed once you have the final results, but it’s wise to keep them until you’ve verified the results.fasta. The per-cluster FASTA files are mainly for debugging or traceability and can also be hefty for large datasets.


---








































