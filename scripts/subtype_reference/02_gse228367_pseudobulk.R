#!/usr/bin/env Rscript

# Biological-replicate count-pseudobulk DCT2-versus-DCT1 analysis.
#
# The author-separated GSE228367 RDS objects are used only for DCT1/DCT2 cell
# membership and biological-replicate metadata. Their RNA@counts slot is
# corrected/fractional and is explicitly rejected as an edgeR input. Cell
# barcodes are instead matched exactly to the official NK1--NK3 filtered 10x
# H5 matrices in GSE228367_RAW.tar. No flight or phosphoproteomic result is
# read.

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
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
  stop("Usage: 02_gse228367_pseudobulk.R --config CONFIG --output-dir DIR")
}
config_path <- args$config
out_dir <- args[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cfg <- yaml.load_file(config_path)
if (!isTRUE(cfg$flight_blind)) {
  stop("Reference config must declare flight_blind: true")
}
rb <- cfg$reference_builder

input_for_role <- function(dataset, role) {
  entries <- rb$inputs[[dataset]]
  matches <- Filter(function(x) identical(x$role, role), entries)
  if (length(matches) != 1L) {
    stop(sprintf("Expected exactly one %s/%s input in config", dataset, role))
  }
  matches[[1]]
}

load_double_gz_rds <- function(path) {
  con <- gzcon(gzfile(path, "rb"))
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

get_component <- function(object, name) {
  value <- attr(object, name)
  if (!is.null(value)) return(value)
  if (isS4(object) && name %in% slotNames(object)) return(slot(object, name))
  NULL
}

metadata_only <- function(object, subtype, sample_column, diet_column, allowed_diets) {
  metadata <- get_component(object, "meta.data")
  if (is.null(metadata)) stop(subtype, " meta.data not found")
  metadata <- as.data.frame(metadata)
  needed <- c(sample_column, diet_column)
  if (!all(needed %in% names(metadata))) {
    stop(subtype, " metadata lacks: ", paste(setdiff(needed, names(metadata)), collapse = ", "))
  }
  observed_diets <- unique(as.character(metadata[[diet_column]]))
  if (!all(observed_diets %in% allowed_diets)) {
    stop(subtype, " contains non-control diet: ", paste(observed_diets, collapse = ", "))
  }
  cell_id <- rownames(metadata)
  if (is.null(cell_id) || any(!nzchar(cell_id)) || anyDuplicated(cell_id)) {
    stop(subtype, " metadata requires unique non-empty cell identifiers")
  }
  replicate <- as.character(metadata[[sample_column]])
  expected_prefix <- paste0(tolower(replicate), "_")
  prefix_ok <- startsWith(tolower(cell_id), expected_prefix)
  if (!all(prefix_ok)) {
    stop(
      subtype, " cell identifiers do not carry the expected replicate prefix; examples: ",
      paste(head(cell_id[!prefix_ok]), collapse = ", ")
    )
  }
  raw_barcode <- substring(cell_id, nchar(expected_prefix) + 1L)
  if (any(!nzchar(raw_barcode))) stop(subtype, " produced an empty raw barcode")
  data.frame(
    cell_id = cell_id,
    raw_barcode = raw_barcode,
    replicate = replicate,
    subtype = subtype,
    stringsAsFactors = FALSE
  )
}

check_fractional_rna_slot <- function(object, subtype) {
  assays <- get_component(object, "assays")
  if (is.null(assays) || !"RNA" %in% names(assays)) {
    return(list(present = FALSE, fractional_fraction = NA_real_))
  }
  counts <- get_component(assays[["RNA"]], "counts")
  if (is.null(counts)) return(list(present = FALSE, fractional_fraction = NA_real_))
  counts <- as(counts, "dgCMatrix")
  values <- counts@x
  fraction <- if (length(values)) mean(abs(values - round(values)) > 1e-8) else NA_real_
  list(
    present = TRUE,
    n_features = nrow(counts),
    n_cells = ncol(counts),
    fractional_fraction = fraction,
    accepted_for_edger = FALSE
  )
}

dct1_entry <- input_for_role("GSE228367", "control_DCT1_seurat_rds")
dct2_entry <- input_for_role("GSE228367", "control_DCT2_seurat_rds")
raw_entry <- input_for_role("GSE228367", "six_sample_raw_10x_archive")
input_paths <- c(dct1_entry$path, dct2_entry$path, raw_entry$path)
missing <- input_paths[!file.exists(input_paths)]
if (length(missing)) {
  stop("Required GSE228367 reference input(s) missing: ", paste(missing, collapse = ", "))
}

sample_column <- rb$gse228367$replicate_column
diet_column <- rb$gse228367$diet_column
allowed_diets <- as.character(rb$gse228367$allowed_control_values)
min_cells <- as.integer(rb$thresholds$min_cells_per_pseudobulk)
minimum_pairs <- as.integer(rb$thresholds$minimum_biological_replicates)

message("Loading GSE228367 author-separated DCT1 membership object")
dct1_object <- load_double_gz_rds(dct1_entry$path)
message("Loading GSE228367 author-separated DCT2 membership object")
dct2_object <- load_double_gz_rds(dct2_entry$path)
dct1_membership <- metadata_only(
  dct1_object, "DCT1", sample_column, diet_column, allowed_diets
)
dct2_membership <- metadata_only(
  dct2_object, "DCT2", sample_column, diet_column, allowed_diets
)
slot_qc <- list(
  DCT1_RNA_counts = check_fractional_rna_slot(dct1_object, "DCT1"),
  DCT2_RNA_counts = check_fractional_rna_slot(dct2_object, "DCT2")
)
rm(dct1_object, dct2_object)
gc()

membership <- rbind(dct1_membership, dct2_membership)
if (anyDuplicated(membership[, c("replicate", "raw_barcode")])) {
  stop("A raw barcode is assigned to both subtype membership objects")
}
membership_path <- file.path(out_dir, "gse228367_subtype_cell_membership.tsv.gz")
write.table(
  membership,
  gzfile(membership_path),
  sep = "\t", quote = FALSE, row.names = FALSE
)

extractor <- rb$gse228367$raw_count_extractor
python <- rb$gse228367$python_executable
if (is.null(extractor) || !file.exists(extractor)) {
  stop("Configured raw-count extractor is missing: ", extractor)
}
if (is.null(python) || !file.exists(python)) {
  stop("Configured Python executable is missing: ", python)
}
extract_args <- c(
  extractor,
  "--archive", raw_entry$path,
  "--membership", membership_path,
  "--output-dir", out_dir
)
status <- system2(python, args = extract_args)
if (!identical(status, 0L)) {
  stop("Official 10x barcode-matching helper failed with status ", status)
}

raw_qc_path <- file.path(out_dir, "gse228367_raw_mapping_qc.json")
raw_qc <- fromJSON(raw_qc_path, simplifyVector = TRUE)
if (!isTRUE(raw_qc$all_barcodes_matched) ||
    !isTRUE(raw_qc$all_counts_integer_like) ||
    !isTRUE(all.equal(raw_qc$minimum_mapping_fraction, 1))) {
  stop("Raw-count mapping QC did not establish exact integer-count reconstruction")
}

counts_table <- read.delim(
  gzfile(file.path(out_dir, "gse228367_raw_pseudobulk_counts.tsv.gz")),
  check.names = FALSE
)
detection_table <- read.delim(
  gzfile(file.path(out_dir, "gse228367_raw_pseudobulk_detection.tsv.gz")),
  check.names = FALSE
)
if (!identical(counts_table$gene_symbol, detection_table$gene_symbol)) {
  stop("Raw pseudobulk count and detection gene orders differ")
}
genes <- as.character(counts_table$gene_symbol)
counts <- as.matrix(counts_table[, -1L, drop = FALSE])
detection <- as.matrix(detection_table[, -1L, drop = FALSE])
rownames(counts) <- genes
rownames(detection) <- genes
storage.mode(counts) <- "numeric"
storage.mode(detection) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0) ||
    max(abs(counts - round(counts))) > 1e-8) {
  stop("Reconstructed pseudobulk matrix is not finite non-negative integer counts")
}
counts <- round(counts)

replicates_dct1 <- sub("^DCT1::", "", grep("^DCT1::", colnames(counts), value = TRUE))
replicates_dct2 <- sub("^DCT2::", "", grep("^DCT2::", colnames(counts), value = TRUE))
paired_reps <- intersect(replicates_dct1, replicates_dct2)
if (length(paired_reps) < minimum_pairs) {
  stop("Insufficient paired DCT1/DCT2 biological replicates")
}
ordered_columns <- c(paste0("DCT1::", paired_reps), paste0("DCT2::", paired_reps))
counts <- counts[, ordered_columns, drop = FALSE]
detection <- detection[, ordered_columns, drop = FALSE]
mapping_qc <- read.delim(
  file.path(out_dir, "gse228367_barcode_mapping_qc.tsv"),
  check.names = FALSE
)
mapping_qc <- mapping_qc[match(ordered_columns, mapping_qc$sample_id), , drop = FALSE]
if (anyNA(mapping_qc$sample_id) || any(mapping_qc$n_cells_membership < min_cells)) {
  stop("Missing or undersized reconstructed pseudobulk")
}

metadata <- data.frame(
  sample_id = ordered_columns,
  replicate = factor(rep(paired_reps, 2L), levels = paired_reps),
  subtype = factor(
    rep(c("DCT1", "DCT2"), each = length(paired_reps)),
    levels = c("DCT1", "DCT2")
  ),
  n_cells = mapping_qc$n_cells_membership,
  barcode_mapping_fraction = mapping_qc$mapping_fraction,
  library_size = colSums(counts),
  count_source = "official_filtered_10x_H5",
  stringsAsFactors = FALSE
)
write.table(
  metadata,
  file.path(out_dir, "gse228367_pseudobulk_sample_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(gene_symbol = rownames(counts), counts, check.names = FALSE),
  gzfile(file.path(out_dir, "gse228367_pseudobulk_counts.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

design <- model.matrix(~ replicate + subtype, metadata)
y <- DGEList(counts = counts)
keep <- filterByExpr(y, design = design)
if (sum(keep) < 100L) {
  stop("Fewer than 100 genes passed pseudobulk expression filtering")
}
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
y <- estimateDisp(y, design = design, robust = TRUE)
fit <- glmQLFit(y, design = design, robust = TRUE)
coef_name <- "subtypeDCT2"
if (!coef_name %in% colnames(design)) {
  stop("Paired design did not produce subtypeDCT2 coefficient")
}
test <- glmQLFTest(fit, coef = which(colnames(design) == coef_name))
de <- topTags(test, n = Inf, sort.by = "none")$table
de$gene_symbol <- rownames(de)

log_cpm <- cpm(y, log = TRUE, prior.count = 0.5)
pair_diffs <- sapply(paired_reps, function(rep_name) {
  log_cpm[, paste0("DCT2::", rep_name)] - log_cpm[, paste0("DCT1::", rep_name)]
})
if (is.null(dim(pair_diffs))) pair_diffs <- matrix(pair_diffs, ncol = 1L)
colnames(pair_diffs) <- paired_reps
n_positive <- rowSums(pair_diffs > 0)
n_negative <- rowSums(pair_diffs < 0)
de$n_consistent_pairs <- ifelse(de$logFC >= 0, n_positive, n_negative)
de$n_pairs <- length(paired_reps)

de$pct_detected_dct1 <- rowMeans(
  detection[de$gene_symbol, paste0("DCT1::", paired_reps), drop = FALSE]
)
de$pct_detected_dct2 <- rowMeans(
  detection[de$gene_symbol, paste0("DCT2::", paired_reps), drop = FALSE]
)
de$log2_fc_dct2_vs_dct1 <- de$logFC
de$fdr <- de$FDR
de <- de[, c(
  "gene_symbol", "log2_fc_dct2_vs_dct1", "logCPM", "F", "PValue", "fdr",
  "pct_detected_dct1", "pct_detected_dct2", "n_consistent_pairs", "n_pairs"
)]
de <- de[order(de$fdr, -abs(de$log2_fc_dct2_vs_dct1)), ]
write.table(
  de,
  file.path(out_dir, "gse228367_dct2_vs_dct1_paired_pseudobulk.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(gene_symbol = rownames(pair_diffs), pair_diffs, check.names = FALSE),
  gzfile(file.path(out_dir, "gse228367_pairwise_logcpm_differences.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

markers <- unique(unlist(rb$qc_marker_panels, use.names = FALSE))
write.table(
  de[de$gene_symbol %in% markers, ],
  file.path(out_dir, "gse228367_marker_qc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

qc <- list(
  dataset = "GSE228367",
  estimand = "DCT2 minus DCT1 in normal-potassium control nuclei",
  model = "~ biological_replicate + subtype; edgeR quasi-likelihood; TMM normalization",
  count_source = "official NK1-NK3 filtered 10x H5 matrices",
  membership_source = "author-separated DCT1/DCT2 RDS metadata",
  rds_RNA_counts_slot_qc = slot_qc,
  rds_counts_used_for_edger = FALSE,
  all_barcodes_matched = raw_qc$all_barcodes_matched,
  minimum_barcode_mapping_fraction = raw_qc$minimum_mapping_fraction,
  all_counts_integer_like = raw_qc$all_counts_integer_like,
  n_paired_biological_replicates = length(paired_reps),
  paired_replicates = paired_reps,
  n_gene_symbols = nrow(counts),
  n_genes_tested = nrow(de),
  min_cells_per_pseudobulk = min_cells,
  all_inputs_reference_only = TRUE
)
write_json(
  qc,
  file.path(out_dir, "gse228367_pseudobulk_qc.json"),
  pretty = TRUE, auto_unbox = TRUE
)

provenance <- list(
  freeze_id = cfg$freeze_id,
  config = config_path,
  config_sha256 = digest(file = config_path, algo = "sha256", serialize = FALSE),
  inputs = list(
    list(path = dct1_entry$path, expected_sha256 = dct1_entry$sha256, role = "membership_only"),
    list(path = dct2_entry$path, expected_sha256 = dct2_entry$sha256, role = "membership_only"),
    list(
      path = raw_entry$path,
      expected_sha256 = raw_entry$sha256,
      actual_sha256 = raw_qc$archive_sha256,
      role = "official_integer_counts"
    )
  ),
  barcode_mapping_qc = raw_qc_path,
  flight_result_inputs_used = list(),
  software = list(
    R = R.version.string,
    edgeR = as.character(packageVersion("edgeR")),
    Matrix = as.character(packageVersion("Matrix"))
  )
)
write_json(
  provenance,
  file.path(out_dir, "gse228367_pseudobulk_provenance.json"),
  pretty = TRUE, auto_unbox = TRUE
)
message("Completed GSE228367 official-count paired pseudobulk analysis: ", out_dir)
