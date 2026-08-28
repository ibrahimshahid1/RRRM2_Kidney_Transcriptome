#!/usr/bin/env Rscript

# Count-pseudobulk validation of DCT1/DCT2/CNT in the author-integrated
# GSE150338 object. Cluster labels are frozen in the config and originate from
# the repository's pre-existing, marker-checked Chen annotation; no flight
# result is read or used.

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
  stop("Usage: 02c_gse150338_fine_subtype_pseudobulk.R --config CONFIG --output-dir DIR")
}
config_path <- args$config
out_dir <- args[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cfg <- yaml.load_file(config_path)
if (!isTRUE(cfg$flight_blind)) stop("Config must declare flight_blind: true")
rb <- cfg$reference_builder
settings <- rb$gse150338

entries <- rb$inputs$GSE150338
matches <- Filter(function(x) identical(x$role, "author_integrated_seurat_rdata"), entries)
if (length(matches) != 1L) stop("Expected one author_integrated_seurat_rdata input")
entry <- matches[[1]]
input_path <- entry$path
if (!file.exists(input_path)) stop("Missing author-integrated object: ", input_path)

message("Loading author-integrated GSE150338 object")
env <- new.env(parent = emptyenv())
con <- gzcon(gzfile(input_path, "rb"))
loaded <- load(con, envir = env)
close(con)
object_name <- settings$integrated_object_name
if (!object_name %in% loaded) stop("Integrated object not found: ", object_name)
object <- env[[object_name]]
rm(env)
gc()

get_component <- function(object, name) {
  value <- attr(object, name)
  if (!is.null(value)) return(value)
  if (isS4(object) && name %in% slotNames(object)) return(slot(object, name))
  NULL
}
assays <- get_component(object, "assays")
metadata <- as.data.frame(get_component(object, "meta.data"))
if (is.null(assays) || !"RNA" %in% names(assays)) stop("RNA assay missing")
counts <- get_component(assays[["RNA"]], "counts")
if (is.null(counts)) stop("Raw RNA counts missing")
counts <- as(counts, "dgCMatrix")
if (!all(c(settings$replicate_column, settings$cluster_column) %in% names(metadata))) {
  stop("Required replicate/cluster metadata fields are absent")
}
if (!is.null(colnames(counts)) && !is.null(rownames(metadata))) {
  positions <- match(colnames(counts), rownames(metadata))
  if (anyNA(positions)) stop("Counts and metadata cell identifiers do not align")
  metadata <- metadata[positions, , drop = FALSE]
}

cluster_to_subtype <- unlist(lapply(names(settings$cluster_map), function(subtype) {
  setNames(rep(subtype, length(settings$cluster_map[[subtype]])), settings$cluster_map[[subtype]])
}))
metadata$subtype <- cluster_to_subtype[as.character(metadata[[settings$cluster_column]])]
keep <- !is.na(metadata$subtype)
counts <- counts[, keep, drop = FALSE]
metadata <- metadata[keep, , drop = FALSE]
metadata$subtype <- factor(metadata$subtype, levels = c("DCT1", "DCT2", "CNT"))
metadata$replicate <- as.character(metadata[[settings$replicate_column]])

min_cells <- as.integer(rb$thresholds$min_cells_per_pseudobulk)
cell_group <- paste(metadata$subtype, metadata$replicate, sep = "::")
levels_group <- as.vector(outer(c("DCT1", "DCT2", "CNT"), sort(unique(metadata$replicate)), paste, sep = "::"))
levels_group <- levels_group[levels_group %in% cell_group]
group_factor <- factor(cell_group, levels = levels_group)
membership <- sparseMatrix(
  i = seq_along(group_factor),
  j = as.integer(group_factor),
  x = 1,
  dims = c(length(group_factor), length(levels_group))
)
pb <- counts %*% membership
colnames(pb) <- levels_group
n_cells <- tabulate(as.integer(group_factor), nbins = length(levels_group))
if (any(n_cells < min_cells)) {
  stop("Fine-subtype pseudobulk below min_cells: ",
       paste(paste(levels_group[n_cells < min_cells], n_cells[n_cells < min_cells], sep = ":"), collapse = ", "))
}
detected <- sweep(as.matrix((counts > 0) %*% membership), 2L, n_cells, "/")
colnames(detected) <- levels_group

sample_parts <- do.call(rbind, strsplit(levels_group, "::", fixed = TRUE))
sample_meta <- data.frame(
  sample_id = levels_group,
  subtype = factor(sample_parts[, 1], levels = c("DCT1", "DCT2", "CNT")),
  replicate = factor(sample_parts[, 2]),
  n_cells = n_cells,
  library_size = Matrix::colSums(pb),
  stringsAsFactors = FALSE
)
complete_reps <- Reduce(intersect, lapply(levels(sample_meta$subtype), function(subtype) {
  as.character(sample_meta$replicate[sample_meta$subtype == subtype])
}))
minimum_reps <- as.integer(rb$thresholds$minimum_biological_replicates)
if (length(complete_reps) < minimum_reps) {
  stop("Fewer than configured complete DCT1/DCT2/CNT replicate blocks")
}
keep_samples <- as.character(sample_meta$replicate) %in% complete_reps
pb <- pb[, keep_samples, drop = FALSE]
detected <- detected[, keep_samples, drop = FALSE]
sample_meta <- sample_meta[keep_samples, , drop = FALSE]
write.table(sample_meta, file.path(out_dir, "gse150338_fine_pseudobulk_sample_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_symbol = rownames(pb), as.matrix(pb), check.names = FALSE),
            gzfile(file.path(out_dir, "gse150338_fine_pseudobulk_counts.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

design <- model.matrix(~ replicate + subtype, sample_meta)
y <- DGEList(pb)
keep_genes <- filterByExpr(y, design = design)
y <- calcNormFactors(y[keep_genes, , keep.lib.sizes = FALSE])
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

coef_dct2 <- which(colnames(design) == "subtypeDCT2")
coef_cnt <- which(colnames(design) == "subtypeCNT")
if (length(coef_dct2) != 1L || length(coef_cnt) != 1L) {
  stop("Expected DCT2 and CNT design coefficients are missing")
}
contrast <- function(dct2 = 0, cnt = 0) {
  v <- numeric(ncol(design))
  v[coef_dct2] <- dct2
  v[coef_cnt] <- cnt
  v
}
run_test <- function(vector) {
  tab <- topTags(glmQLFTest(fit, contrast = vector), n = Inf, sort.by = "none")$table
  tab$gene_symbol <- rownames(tab)
  tab
}
t_dct2_dct1 <- run_test(contrast(dct2 = 1))
t_dct2_cnt <- run_test(contrast(dct2 = 1, cnt = -1))
t_dct1_cnt <- run_test(contrast(cnt = -1))

genes <- t_dct2_dct1$gene_symbol
log_cpm <- cpm(y, log = TRUE, prior.count = 0.5)
pair_diffs <- sapply(complete_reps, function(rep_name) {
  log_cpm[, paste0("DCT2::", rep_name)] - log_cpm[, paste0("DCT1::", rep_name)]
})
if (is.null(dim(pair_diffs))) pair_diffs <- matrix(pair_diffs, ncol = 1L)
n_consistent <- ifelse(
  t_dct2_dct1$logFC >= 0,
  rowSums(pair_diffs > 0),
  rowSums(pair_diffs < 0)
)
detection_for <- function(subtype) {
  cols <- paste0(subtype, "::", complete_reps)
  rowMeans(detected[rownames(y), cols, drop = FALSE])
}
result <- data.frame(
  gene_symbol = genes,
  log2_fc_dct2_vs_dct1 = t_dct2_dct1$logFC,
  fdr_dct2_vs_dct1 = t_dct2_dct1$FDR,
  log2_fc_dct2_vs_cnt = t_dct2_cnt$logFC,
  fdr_dct2_vs_cnt = t_dct2_cnt$FDR,
  log2_fc_dct1_vs_cnt = t_dct1_cnt$logFC,
  fdr_dct1_vs_cnt = t_dct1_cnt$FDR,
  pct_detected_dct1 = detection_for("DCT1"),
  pct_detected_dct2 = detection_for("DCT2"),
  pct_detected_cnt = detection_for("CNT"),
  n_consistent_dct2_vs_dct1 = n_consistent,
  n_pairs = length(complete_reps),
  stringsAsFactors = FALSE
)
result <- result[order(result$fdr_dct2_vs_dct1, -abs(result$log2_fc_dct2_vs_dct1)), ]
write.table(result, file.path(out_dir, "gse150338_fine_subtype_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

marker_genes <- unique(unlist(rb$qc_marker_panels, use.names = FALSE))
marker_qc <- result[result$gene_symbol %in% marker_genes, ]
write.table(marker_qc, file.path(out_dir, "gse150338_fine_marker_qc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

qc <- list(
  dataset = "GSE150338",
  source_url = entry$source_url,
  input_sha256_expected = entry$sha256,
  input_sha256_actual = digest(file = input_path, algo = "sha256", serialize = FALSE),
  integrated_object = object_name,
  raw_count_assay = "RNA",
  replicate_column = settings$replicate_column,
  cluster_column = settings$cluster_column,
  cluster_map = settings$cluster_map,
  cluster_map_source = settings$cluster_map_source,
  n_cells_retained = ncol(counts),
  n_complete_replicates = length(complete_reps),
  complete_replicates = complete_reps,
  cell_counts = as.list(table(metadata$subtype, metadata$replicate)),
  model = "~ biological_replicate + subtype; edgeR quasi-likelihood; TMM normalization",
  flight_result_inputs_used = list()
)
write_json(qc, file.path(out_dir, "gse150338_fine_subtype_provenance.json"),
           pretty = TRUE, auto_unbox = TRUE)
message("Completed GSE150338 fine-subtype count-pseudobulk validation: ", out_dir)
