#!/usr/bin/env Rscript
# Multi-study residualization + harmonization (Phase 1 of the Cross-OSDR framework).

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(limma)
})

option_list <- list(
    make_option("--counts", type = "character", help = "Path to raw or VST counts (genes x samples)."),
    make_option("--meta", type = "character", help = "Path to per-study metadata table."),
    make_option("--clr", type = "character", default = "",
                help = "Optional CLR-transformed cell composition table (samples x cell types)."),
    make_option("--gene_universe", type = "character", default = "",
                help = "Optional common gene universe TSV; restricts output to these Ensembl IDs."),
    make_option("--study", type = "character", help = "Study ID, e.g. OSD-771, OSD-513, OSD-253."),
    make_option("--outdir", type = "character", help = "Output directory for residuals + manifest."),
    make_option("--tech_cols", type = "character", default = "LibraryBatch,ReadDepth,RIN",
                help = "Comma-separated technical covariate column names in metadata."),
    make_option("--sample_col", type = "character", default = "Sample Name",
                help = "Sample identifier column name in metadata."),
    make_option("--force_keep", type = "character", default = "",
                help = "Optional file with Ensembl IDs that must be kept regardless of gene_universe.")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$counts), !is.null(opt$meta), !is.null(opt$study), !is.null(opt$outdir))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

strip_ensembl_version <- function(x) sub("\\.\\d+$", "", as.character(x))

cat(sprintf("[%s] Loading counts: %s\n", opt$study, opt$counts))
counts_dt <- fread(opt$counts)
gene_col <- names(counts_dt)[1]
gene_ids <- strip_ensembl_version(counts_dt[[gene_col]])
counts <- as.matrix(counts_dt[, -1, with = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- gene_ids
if (anyDuplicated(rownames(counts)) > 0) {
    counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
}

cat(sprintf("[%s] Loading metadata: %s\n", opt$study, opt$meta))
meta <- as.data.frame(fread(opt$meta))
if (!opt$sample_col %in% colnames(meta)) {
    stop(sprintf("Sample column %s not in metadata", opt$sample_col))
}
common <- intersect(colnames(counts), meta[[opt$sample_col]])
if (length(common) == 0) {
    stop("No overlapping samples between counts and metadata")
}
counts <- counts[, common, drop = FALSE]
meta <- meta[match(common, meta[[opt$sample_col]]), , drop = FALSE]
cat(sprintf("[%s] Samples retained: %d\n", opt$study, length(common)))

design_terms <- c()
covariates <- list()

if (nzchar(opt$clr)) {
    cat(sprintf("[%s] Loading cell-composition CLR: %s\n", opt$study, opt$clr))
    clr <- as.data.frame(fread(opt$clr))
    clr_sample_col <- names(clr)[1]
    rownames(clr) <- clr[[clr_sample_col]]
    clr <- clr[common, , drop = FALSE]
    clr_terms <- setdiff(colnames(clr), clr_sample_col)
    for (term in clr_terms) {
        covariates[[term]] <- as.numeric(clr[[term]])
        design_terms <- c(design_terms, term)
    }
}

tech_cols <- strsplit(opt$tech_cols, ",", fixed = TRUE)[[1]]
for (col in tech_cols) {
    if (col %in% colnames(meta)) {
        val <- meta[[col]]
        if (is.character(val) || is.factor(val)) {
            covariates[[col]] <- factor(val)
        } else {
            covariates[[col]] <- as.numeric(val)
        }
        design_terms <- c(design_terms, col)
    } else {
        cat(sprintf("[%s] Skipping missing tech covariate: %s\n", opt$study, col))
    }
}

if (length(design_terms) == 0) {
    stop(sprintf("[%s] No covariates resolved; refusing to write null residualization.", opt$study))
}

design_df <- as.data.frame(covariates)
design_df[is.na(design_df)] <- 0
design_formula <- as.formula(paste("~", paste(design_terms, collapse = "+")))
design <- model.matrix(design_formula, data = design_df)
cat(sprintf("[%s] Design matrix: %d samples x %d columns (terms=%s)\n",
            opt$study, nrow(design), ncol(design), paste(design_terms, collapse = ", ")))

if (max(counts) > 30) {
    cat(sprintf("[%s] Counts look like raw integers; applying CPM + log2 before residualization.\n", opt$study))
    counts_norm <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
} else {
    cat(sprintf("[%s] Counts look already normalized; using as-is.\n", opt$study))
    counts_norm <- counts
}

fit <- limma::lmFit(counts_norm, design = design)
residuals_mat <- counts_norm - fit$coefficients %*% t(design)

force_keep_ids <- c()
if (nzchar(opt$force_keep) && file.exists(opt$force_keep)) {
    force_keep_ids <- strip_ensembl_version(readLines(opt$force_keep))
}

if (nzchar(opt$gene_universe) && file.exists(opt$gene_universe)) {
    keep_ids <- strip_ensembl_version(readLines(opt$gene_universe))
    if (length(force_keep_ids) > 0) keep_ids <- union(keep_ids, force_keep_ids)
    keep_ids <- intersect(keep_ids, rownames(residuals_mat))
    residuals_mat <- residuals_mat[keep_ids, , drop = FALSE]
    cat(sprintf("[%s] Restricted to common gene universe: %d genes\n", opt$study, nrow(residuals_mat)))
}

residuals_path <- file.path(opt$outdir, "residuals.tsv.gz")
con <- gzfile(residuals_path, "w")
write.table(
    data.frame(gene_id = rownames(residuals_mat), residuals_mat, check.names = FALSE),
    con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
close(con)
cat(sprintf("[%s] Wrote residuals: %s\n", opt$study, residuals_path))

manifest <- list(
    study = opt$study,
    counts = opt$counts,
    metadata = opt$meta,
    clr = opt$clr,
    gene_universe = opt$gene_universe,
    formula = deparse(design_formula),
    design_terms = design_terms,
    n_samples = length(common),
    n_genes_in = nrow(counts_norm),
    n_genes_out = nrow(residuals_mat),
    residuals = residuals_path,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)
writeLines(
    jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE),
    con = file.path(opt$outdir, "manifest.json")
)
cat(sprintf("[%s] Manifest: %s\n", opt$study, file.path(opt$outdir, "manifest.json")))
