#!/usr/bin/env Rscript

# MuSiC reference-atlas sensitivity analysis.
#
# This script compares primary MuSiC cell-type proportions against an alternate
# deconvolution/reference-free estimate (for example bMIND, CIBERSORTx-compatible
# output, or a TMS pseudo-bulk holdout). It also optionally compares key
# downstream Phase 3 top-decile rewiring and Phase 7 pathway summaries produced
# after rerunning residualization with alternate proportions.

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    music = "data/processed/deconvolution/music_segment_direct_proportions.csv",
    alternate = "",
    primary_phase3 = "",
    alternate_phase3 = "",
    primary_phase7 = "",
    alternate_phase7 = "",
    outdir = "data/results/deconvolution_sensitivity",
    top_quantile = 0.90
  )
  for (arg in args) {
    if (grepl("^--", arg)) {
      parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- if (length(parts) > 1) parts[[2]] else TRUE
      out[[key]] <- value
    }
  }
  out$top_quantile <- as.numeric(out$top_quantile)
  out
}

read_table_auto <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(sprintf("Missing table: %s", path))
  }
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  read.table(path, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

sample_col <- function(df) {
  hits <- intersect(c("sample", "Sample", "Sample Name", "Sample Name (raw_counts_colname)"), names(df))
  if (length(hits) > 0) hits[[1]] else names(df)[[1]]
}

compare_proportions <- function(primary, alternate) {
  s1 <- sample_col(primary)
  s2 <- sample_col(alternate)
  names(primary)[names(primary) == s1] <- "sample"
  names(alternate)[names(alternate) == s2] <- "sample"
  merged <- merge(primary, alternate, by = "sample", suffixes = c("_music", "_alternate"))
  shared <- intersect(sub("_music$", "", grep("_music$", names(merged), value = TRUE)),
                      sub("_alternate$", "", grep("_alternate$", names(merged), value = TRUE)))
  rows <- list()
  for (cell in shared) {
    x <- merged[[paste0(cell, "_music")]]
    y <- merged[[paste0(cell, "_alternate")]]
    rows[[length(rows) + 1]] <- data.frame(
      cell_type = cell,
      n_samples = length(x),
      pearson = suppressWarnings(cor(x, y, method = "pearson", use = "pairwise.complete.obs")),
      spearman = suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs")),
      mean_abs_delta = mean(abs(x - y), na.rm = TRUE)
    )
  }
  do.call(rbind, rows)
}

top_decile_jaccard <- function(primary_path, alternate_path, q = 0.90) {
  p <- read_table_auto(primary_path)
  a <- read_table_auto(alternate_path)
  stat_col <- function(df) {
    hits <- intersect(c("rewiring_mean", "edge_sum_node_rewiring_obs", "rewiring"), names(df))
    if (length(hits) == 0) stop("No rewiring statistic column found")
    hits[[1]]
  }
  ps <- stat_col(p)
  as <- stat_col(a)
  p_top <- p$gene[p[[ps]] >= quantile(p[[ps]], q, na.rm = TRUE)]
  a_top <- a$gene[a[[as]] >= quantile(a[[as]], q, na.rm = TRUE)]
  data.frame(
    primary_file = basename(primary_path),
    alternate_file = basename(alternate_path),
    top_quantile = q,
    n_primary = length(p_top),
    n_alternate = length(a_top),
    jaccard = length(intersect(p_top, a_top)) / length(union(p_top, a_top))
  )
}

pathway_overlap <- function(primary_path, alternate_path) {
  p <- read_table_auto(primary_path)
  a <- read_table_auto(alternate_path)
  feature_col <- intersect(c("pathway", "gene_set", "term"), names(p))
  if (length(feature_col) == 0) stop("Primary pathway table lacks pathway/gene_set/term column")
  col <- feature_col[[1]]
  if (!col %in% names(a)) stop("Alternate pathway table lacks matching feature column")
  p_sig <- p[[col]]
  a_sig <- a[[col]]
  data.frame(
    primary_file = basename(primary_path),
    alternate_file = basename(alternate_path),
    n_primary = length(unique(p_sig)),
    n_alternate = length(unique(a_sig)),
    jaccard = length(intersect(unique(p_sig), unique(a_sig))) / length(union(unique(p_sig), unique(a_sig)))
  )
}

main <- function() {
  args <- parse_args()
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  if (!nzchar(args$alternate)) {
    stop("Provide --alternate=<alternate proportion table>. This may be bMIND, CIBERSORTx-compatible, or pseudo-bulk holdout output.")
  }

  music <- read_table_auto(args$music)
  alt <- read_table_auto(args$alternate)
  prop_summary <- compare_proportions(music, alt)
  write.table(prop_summary, file = file.path(args$outdir, "proportion_sensitivity.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  stability_rows <- list()
  if (nzchar(args$primary_phase3) && nzchar(args$alternate_phase3)) {
    stability_rows[[length(stability_rows) + 1]] <- top_decile_jaccard(
      args$primary_phase3, args$alternate_phase3, args$top_quantile
    )
  }
  if (length(stability_rows) > 0) {
    write.table(do.call(rbind, stability_rows), file = file.path(args$outdir, "top_decile_jaccard.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  if (nzchar(args$primary_phase7) && nzchar(args$alternate_phase7)) {
    path_summary <- pathway_overlap(args$primary_phase7, args$alternate_phase7)
    write.table(path_summary, file = file.path(args$outdir, "pathway_stability.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  cat("[OK] MuSiC deconvolution sensitivity outputs written to", args$outdir, "\n")
}

main()
