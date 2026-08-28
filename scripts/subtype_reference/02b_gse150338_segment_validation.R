#!/usr/bin/env Rscript

# Independent segment-level DCT-to-CNT retention validation from GSE150338.
# The source is processed replicate TPM, so this is not a count-DE analysis.
# Technical a/b/c aliquots are averaged within mouse before paired inference.

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(digest)
})

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--") || i == length(x)) {
      stop("Arguments must be supplied as --key value")
    }
    out[[substring(key, 3L)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$config) || is.null(args[["output-dir"]])) {
  stop("Usage: 02b_gse150338_segment_validation.R --config CONFIG --output-dir DIR")
}
config_path <- args$config
out_dir <- args[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cfg <- yaml.load_file(config_path)
if (!isTRUE(cfg$flight_blind)) {
  stop("Reference config must declare flight_blind: true")
}
rb <- cfg$reference_builder

entries <- rb$inputs$GSE150338
matches <- Filter(function(x) identical(x$role, "microdissected_segment_replicate_TPM"), entries)
if (length(matches) != 1L) {
  stop("Config must define exactly one GSE150338 microdissected_segment_replicate_TPM input")
}
entry <- matches[[1]]
input_path <- entry$path
if (!file.exists(input_path)) {
  stop("Missing GSE150338 segment TPM file: ", input_path)
}

settings <- rb$gse150338
margin <- as.numeric(rb$thresholds$external_cnt_retention_log2_margin)
fdr_cutoff <- as.numeric(rb$thresholds$external_retention_fdr)
detection_tpm <- as.numeric(settings$tpm_detection_threshold)
pseudocount <- as.numeric(settings$log2_pseudocount)

raw <- read.delim(
  gzfile(input_path),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)
gene_col <- settings$gene_symbol_column
if (!gene_col %in% names(raw)) {
  stop("Gene-symbol column absent: ", gene_col)
}
pattern <- settings$biological_replicate_regex
replicate_columns <- grep(pattern, names(raw), value = TRUE)
if (!length(replicate_columns)) {
  stop("No DCT/CNT/CCD biological-replicate columns matched configured regex")
}

parsed <- regexec(pattern, replicate_columns)
parts <- regmatches(replicate_columns, parsed)
sample_map <- data.frame(
  column = replicate_columns,
  mouse = vapply(parts, `[[`, character(1), 2L),
  segment = vapply(parts, `[[`, character(1), 3L),
  stringsAsFactors = FALSE
)
sample_map <- sample_map[sample_map$segment %in% settings$segments, ]

values <- as.matrix(raw[, sample_map$column, drop = FALSE])
storage.mode(values) <- "double"
genes <- trimws(as.character(raw[[gene_col]]))
valid <- !is.na(genes) & nzchar(genes) & rowSums(is.finite(values)) > 0
values <- values[valid, , drop = FALSE]
genes <- genes[valid]
values[!is.finite(values)] <- 0
if (anyDuplicated(genes)) {
  values <- rowsum(values, group = genes, reorder = FALSE)
  genes <- rownames(values)
} else {
  rownames(values) <- genes
}

mouse_segment <- paste(sample_map$mouse, sample_map$segment, sep = "::")
collapsed_names <- unique(mouse_segment)
collapsed <- sapply(collapsed_names, function(label) {
  rowMeans(values[, mouse_segment == label, drop = FALSE])
})
if (is.null(dim(collapsed))) {
  collapsed <- matrix(collapsed, ncol = 1L)
}
rownames(collapsed) <- genes
colnames(collapsed) <- collapsed_names

segment_mice <- lapply(settings$segments, function(segment) {
  sub(paste0("::", segment, "$"), "", grep(paste0("::", segment, "$"), colnames(collapsed), value = TRUE))
})
names(segment_mice) <- settings$segments
paired_dct_cnt <- intersect(segment_mice$DCT, segment_mice$CNT)
paired_cnt_ccd <- intersect(segment_mice$CNT, segment_mice$CCD)
if (length(paired_dct_cnt) < 3L) {
  stop("Fewer than three paired DCT/CNT mice in GSE150338")
}

log_tpm <- log2(collapsed + pseudocount)
dct_cnt_diff <- sapply(paired_dct_cnt, function(mouse) {
  log_tpm[, paste0(mouse, "::CNT")] - log_tpm[, paste0(mouse, "::DCT")]
})
cnt_ccd_diff <- sapply(paired_cnt_ccd, function(mouse) {
  log_tpm[, paste0(mouse, "::CNT")] - log_tpm[, paste0(mouse, "::CCD")]
})
if (is.null(dim(dct_cnt_diff))) dct_cnt_diff <- matrix(dct_cnt_diff, ncol = 1L)
if (is.null(dim(cnt_ccd_diff))) cnt_ccd_diff <- matrix(cnt_ccd_diff, ncol = 1L)

paired_summary <- function(diff, noninferiority_margin = NULL) {
  n <- ncol(diff)
  mean_effect <- rowMeans(diff)
  sd_effect <- apply(diff, 1L, sd)
  se <- sd_effect / sqrt(n)
  t_two <- mean_effect / se
  p_two <- 2 * pt(-abs(t_two), df = n - 1L)
  p_two[!is.finite(p_two)] <- ifelse(mean_effect[!is.finite(p_two)] == 0, 1, 0)
  out <- data.frame(
    mean_effect = mean_effect,
    se = se,
    p_two_sided = p_two,
    fdr_two_sided = p.adjust(p_two, method = "BH")
  )
  if (!is.null(noninferiority_margin)) {
    t_ni <- (mean_effect - noninferiority_margin) / se
    p_ni <- pt(t_ni, df = n - 1L, lower.tail = FALSE)
    p_ni[!is.finite(p_ni)] <- ifelse(
      mean_effect[!is.finite(p_ni)] > noninferiority_margin, 0, 1
    )
    out$ci_low_one_sided_95 <- mean_effect - qt(0.95, df = n - 1L) * se
    out$p_noninferiority <- p_ni
    out$fdr_noninferiority <- p.adjust(p_ni, method = "BH")
  }
  out
}

dct_cnt <- paired_summary(dct_cnt_diff, margin)
cnt_ccd <- paired_summary(cnt_ccd_diff)
cnt_columns <- paste0(segment_mice$CNT, "::CNT")
cnt_detect <- rowMeans(collapsed[, cnt_columns, drop = FALSE] >= detection_tpm)

result <- data.frame(
  gene_symbol = rownames(collapsed),
  log2_fc_cnt_vs_dct = dct_cnt$mean_effect,
  ci_low_cnt_vs_dct = dct_cnt$ci_low_one_sided_95,
  p_cnt_vs_dct = dct_cnt$p_two_sided,
  fdr_cnt_vs_dct = dct_cnt$fdr_two_sided,
  p_noninferiority_cnt_vs_dct = dct_cnt$p_noninferiority,
  fdr_noninferiority_cnt_vs_dct = dct_cnt$fdr_noninferiority,
  log2_fc_cnt_vs_ccd = cnt_ccd$mean_effect,
  fdr_cnt_vs_ccd = cnt_ccd$fdr_two_sided,
  pct_replicates_detected_cnt = cnt_detect,
  n_consistent_cnt_retained = rowSums(dct_cnt_diff > margin),
  n_pairs = length(paired_dct_cnt),
  stringsAsFactors = FALSE
)
result <- result[order(result$fdr_noninferiority_cnt_vs_dct, -result$log2_fc_cnt_vs_dct), ]
write.table(
  result,
  file.path(out_dir, "gse150338_segment_transition_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  sample_map,
  file.path(out_dir, "gse150338_segment_sample_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

qc <- list(
  dataset = "GSE150338",
  source_url = entry$source_url,
  input_sha256_expected = entry$sha256,
  input_sha256_actual = digest(file = input_path, algo = "sha256", serialize = FALSE),
  expression_unit = "author-provided TPM",
  technical_aliquots_collapsed_within_mouse = TRUE,
  paired_DCT_CNT_mice = paired_dct_cnt,
  paired_CNT_CCD_mice = paired_cnt_ccd,
  noninferiority_margin_log2_CNT_vs_DCT = margin,
  retention_fdr_cutoff = fdr_cutoff,
  interpretation_limit = paste(
    "Microdissected DCT is not relabeled as DCT1 or DCT2.",
    "This analysis validates DCT-to-CNT retention, not a strict DCT2 peak."
  ),
  flight_result_inputs_used = list()
)
write_json(
  qc,
  file.path(out_dir, "gse150338_segment_validation_provenance.json"),
  pretty = TRUE, auto_unbox = TRUE
)
message("Completed GSE150338 segment validation: ", out_dir)
