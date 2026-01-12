suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(Biobase)
  library(edgeR)
  library(zellkonverter)
  library(MuSiC)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# ---------------------------
# USER PATHS (edit these)
# ---------------------------
bulk_counts_csv <- "data/processed/aligned_outputs/rsem_rRNArm_raw_counts.csv"
bulk_meta_csv   <- "data/processed/aligned_outputs/metadata_aligned.tsv"

# Use the H5AD you wrote (best), NOT the obs/var CSV
# e.g. tms_kidney_female/tms_kidney_female_ALLDATASETS_innerGenes.h5ad
sc_h5ad_path <- "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"


outdir <- "data/processed/deconvolution"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Helpers
# ---------------------------

strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

# ---------------------------
# ID Space Detection & Marker Resolution Helpers
# ---------------------------
detect_id_space <- function(ids) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) return("unknown")
  frac_ens <- mean(grepl("^ENSMUSG", ids))
  if (frac_ens > 0.5) "ensembl" else "symbol"
}

symbols_to_ensembl <- function(symbols) {
  map <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = "ENSEMBL"
  )
  unique(na.omit(map$ENSEMBL))
}

ensembl_to_symbols <- function(ensembl_ids) {
  map <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(ensembl_ids),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  )
  unique(na.omit(map$SYMBOL))
}

resolve_markers_to_matrix <- function(markers, matrix_rownames) {
  space <- detect_id_space(matrix_rownames)
  
  # If user passes Ensembl IDs already
  looks_ens <- grepl("^ENSMUSG", markers)
  
  if (space == "ensembl") {
    ens <- if (all(looks_ens)) unique(markers) else symbols_to_ensembl(markers)
    return(intersect(ens, matrix_rownames))
  } else if (space == "symbol") {
    sym <- if (all(looks_ens)) ensembl_to_symbols(markers) else unique(markers)
    return(intersect(sym, matrix_rownames))
  } else {
    warning("Could not detect ID space of matrix rownames.")
    return(character(0))
  }
}

# Heuristic: decide whether matrix is genes x samples or samples x genes
ensure_genes_by_samples <- function(mat) {
  # If rownames look like sample IDs more than gene IDs, transpose.
  # We'll use a simple heuristic: Ensembl-like rows vs columns.
  r_ens <- mean(grepl("^ENSMUSG", rownames(mat)))
  c_ens <- mean(grepl("^ENSMUSG", colnames(mat)))

  if (!is.na(c_ens) && c_ens > r_ens) {
    message("Transposing bulk matrix (detected genes on columns).")
    mat <- t(mat)
  }
  mat
}

# ---------------------------
# 1) Load bulk counts + metadata
# ---------------------------
bulk_dt <- fread(bulk_counts_csv)
# Assume first column is gene id if it isn't numeric
gene_col <- names(bulk_dt)[1]
bulk_genes <- bulk_dt[[gene_col]]
bulk_mat <- as.matrix(bulk_dt[, -1, with = FALSE])
rownames(bulk_mat) <- bulk_genes

# Check orientation
bulk_mat <- ensure_genes_by_samples(bulk_mat)

# Filter out empty/NA rows
bulk_mat <- bulk_mat[complete.cases(rownames(bulk_mat)), , drop = FALSE]
storage.mode(bulk_mat) <- "numeric"

# Strip Ensembl versions and collapse duplicates
rownames(bulk_mat) <- strip_ensembl_version(rownames(bulk_mat))

if (anyDuplicated(rownames(bulk_mat)) > 0) {
  message("Bulk has duplicate gene IDs after stripping versions. Collapsing by sum...")
  bulk_mat <- rowsum(bulk_mat, group = rownames(bulk_mat), reorder = FALSE)
}
stopifnot(anyDuplicated(rownames(bulk_mat)) == 0)

# --- B. Process Metadata ---
meta <- fread(bulk_meta_csv, header = TRUE)
meta <- as.data.frame(meta)

# Simple metadata handling: use Sample Name directly
stopifnot("Sample Name" %in% names(meta))
stopifnot(identical(meta[["Sample Name"]], colnames(bulk_mat)))
rownames(meta) <- meta[["Sample Name"]]
stopifnot(anyDuplicated(rownames(meta)) == 0)

message("Successfully aligned ", ncol(bulk_mat), " samples.")

# ---------------------------
# 2) Load single-cell reference from H5AD
# ---------------------------
# ---------------------------
# 2) Load single-cell reference from Raw MTX
# ---------------------------
message("Loading single-cell reference from Raw MTX...")
raw_dir <- "data/processed/deconvolution/raw_ref"
mtx_path <- file.path(raw_dir, "matrix.mtx")
obs_path <- file.path(raw_dir, "barcodes.tsv")
var_path <- file.path(raw_dir, "features.tsv")

if (!file.exists(mtx_path)) stop("Raw counts export missing. Run scripts/export_raw_counts.py")

raw_mat <- readMM(mtx_path)
# data.table::fread is more robust for csv/tsv
obs <- data.frame(fread(obs_path, header=TRUE), row.names=1, check.names=FALSE)
var <- data.frame(fread(var_path, header=TRUE), row.names=1, check.names=FALSE)

# Handle potential transpose issues with mtx (Scanpy writes shape n_obs x n_vars, readMM reads as such)
# But SCE expects genes x cells (n_vars x n_obs). Scanpy usually writes cells x genes?
# Creating scipy.io.mmwrite(..., raw_adata.X) usually writes standard MM.
# If raw_adata.X is (cells, genes), then readMM returns (rows=cells, cols=genes).
# SCE needs (rows=genes, cols=cells). So we likely need to Transpose.
# Let's check dimensions carefully.
message("MTX dims: ", paste(dim(raw_mat), collapse=" x "))
message("Obs rows: ", nrow(obs))
message("Var rows: ", nrow(var))

if (ncol(raw_mat) == nrow(var)) {
  message("Transposing MTX to match Genes x Cells...")
  raw_mat <- t(raw_mat)
}

# Verify dimensions
if (nrow(raw_mat) != nrow(var)) stop("Matrix rows != n_genes")
if (ncol(raw_mat) != nrow(obs)) stop("Matrix cols != n_cells")

colnames(raw_mat) <- rownames(obs)
rownames(raw_mat) <- rownames(var)

sce <- SingleCellExperiment(assays = list(counts = raw_mat), colData = obs)

# Ensure we have "counts" (created above)
message("Assays in SCE: ", paste(assayNames(sce), collapse = ", "))



# ---------------------------
# 2b) Choose best label column automatically based on podo/glom content
# ---------------------------
# Find available label columns
label_cols <- intersect(c("free_annotation", "cell_type", "subtissue"), colnames(colData(sce)))
message("Available label columns: ", paste(label_cols, collapse = ", "))

# Score each column based on how many podo/glom labels it contains
score_col <- function(col) {
  x <- as.character(colData(sce)[[col]])
  x <- x[!is.na(x)]
  sum(grepl("podo|glom", x, ignore.case = TRUE))
}

scores <- sapply(label_cols, score_col)
message("\nLabel column scores (podo/glom matches):")
print(scores)

# Select the column with highest score
# But prefer free_annotation if available because our segment_map is built on it
if ("free_annotation" %in% names(scores)) {
  best <- "free_annotation"
} else {
  best <- names(which.max(scores))
}
message("\nUsing label column: ", best)

# Build cell_type_use from the best column
ct_use <- as.character(colData(sce)[[best]])
ct_use[is.na(ct_use) | ct_use == ""] <- "Unknown"
colData(sce)$cell_type_use <- factor(ct_use)

# Show top cell types and Unknown fraction
message("\n--- Cell Type Distribution ---")
tab <- sort(table(colData(sce)$cell_type_use), decreasing = TRUE)
print(head(tab, 30))
cat("\nUnknown fraction:", tab["Unknown"] / sum(tab), "\n")

# Pattern-based checks for podocyte and glomerular cells
cat("Podocyte-like cells (podo pattern):", sum(grepl("podo", ct_use, ignore.case = TRUE)), "\n")
cat("Glomerular cells (glom pattern):", sum(grepl("glom", ct_use, ignore.case = TRUE)), "\n")

message("\nCell types matching podo|glom pattern:")
print(table(colData(sce)$cell_type_use[grepl("podo|glom", colData(sce)$cell_type_use, ignore.case = TRUE)]))

# ----------------------------
# 0) Clean bad labels (drop "nan" etc.)
# ----------------------------
ct <- as.character(colData(sce)$cell_type_use)
ct <- trimws(ct)
ct[tolower(ct) %in% c("", "nan", "na", "null", "unknown")] <- NA

keep <- !is.na(ct)
message("Dropping ", sum(!keep), " cells with bad/NA labels (including 'nan')")
sce <- sce[, keep]
colData(sce)$cell_type_use <- factor(ct[keep])

message("\n--- Cell Type Distribution (after cleaning) ---")
print(sort(table(colData(sce)$cell_type_use), decreasing = TRUE))

# ----------------------------
# 1) Segment mapping (keep what you already implemented)
# ----------------------------
segment_from_ct <- function(x) {
  if (is.na(x)) return(NA_character_)
  x2 <- tolower(x)
  
  if (grepl("podocyte", x2)) return("Podocyte")
  if (grepl("mesangial", x2)) return("Mesangial")
  if (grepl("pecam|endothelial|capillary|artery", x2)) return("Endothelial")
  if (grepl("fibroblast|stroma", x2)) return("Fibroblast")
  if (grepl("^cd45|macrophage|\\bt cell\\b|\\bb cell\\b|plasma cell|nk cell|lymph", x2)) return("Immune")
  if (grepl("proximal", x2)) return("PT")
  if (grepl("distal convoluted", x2)) return("DCT")
  if (grepl("thick ascending|loop of henle|tal", x2)) return("TAL_LOH")
  if (grepl("collecting duct|principal", x2)) return("CD")
  if (grepl("brush cell|tuft cell", x2)) return("Other")
  
  "Other"
}

seg <- vapply(as.character(colData(sce)$cell_type_use), segment_from_ct, character(1))
colData(sce)$segment_use <- factor(seg)

# Remove any NA segment cells and clean up factor levels
sce <- sce[, !is.na(colData(sce)$segment_use)]
colData(sce)$segment_use <- droplevels(colData(sce)$segment_use)

message("\n--- Segment Distribution ---")
print(table(colData(sce)$segment_use, useNA = "ifany"))

# ----------------------------
# 2) Coarse groups for tree-guided MuSiC
# ----------------------------
seg <- as.character(colData(sce)$segment_use)

colData(sce)$clusterType <- ifelse(seg %in% c("Podocyte","Endothelial","Mesangial"), "Glomerular",
                           ifelse(seg %in% c("Immune"), "Immune",
                           ifelse(seg %in% c("Fibroblast"), "Stroma", "Tubule")))
colData(sce)$clusterType <- factor(colData(sce)$clusterType)

clusters.type <- list(
  Glomerular = c("Podocyte","Endothelial","Mesangial"),
  Tubule     = c("PT","TAL_LOH","DCT","CD","Other"),
  Immune     = c("Immune"),
  Stroma     = c("Fibroblast")
)

message("\n--- Coarse Group Distribution ---")
print(table(colData(sce)$clusterType))

# Deduplication function for SCE (sparse-safe)
dedup_sce_genes_sum <- function(sce_obj) {
  g <- sub("\\.\\d+$", "", rownames(sce_obj))
  
  if (anyDuplicated(g) == 0) {
    rownames(sce_obj) <- g
    return(sce_obj)
  }
  
  message("SCE has duplicate gene IDs after stripping versions. Collapsing by sum...")
  
  f <- factor(g, levels = unique(g))
  grp <- as.integer(f)
  
  # Indicator matrix: (n_unique_genes x n_genes)
  G <- sparseMatrix(
    i = grp,
    j = seq_along(grp),
    x = 1,
    dims = c(nlevels(f), length(grp))
  )
  
  new_counts <- G %*% counts(sce_obj)
  
  out <- SingleCellExperiment(
    assays = list(counts = new_counts),
    colData = colData(sce_obj)
  )
  rownames(out) <- levels(f)
  
  out
}

# Apply deduplication
sce <- dedup_sce_genes_sum(sce)
stopifnot(anyDuplicated(rownames(sce)) == 0)

# Optional: quick check if it's wildly non-count-like
vals <- as.numeric(counts(sce)[1:min(1000, nrow(sce)), 1:min(5, ncol(sce))])
if (any(vals < 0, na.rm = TRUE)) stop("Reference assay has negative values; not usable as counts.")
if (mean(abs(vals - round(vals)) < 1e-6, na.rm = TRUE) < 0.2) {
  message("Warning: reference assay looks non-integer (may be normalized). MuSiC can still run, but counts are preferable.")
}

# ---------------------------
# 3) Gene intersection after deduplication
# ---------------------------
common_genes <- intersect(rownames(bulk_mat), rownames(sce))
message("Common genes: ", length(common_genes))

if (length(common_genes) < 3000) {
  message("Warning: low overlap. Likely Ensembl vs symbol mismatch. You may need mapping (org.Mm.eg.db / biomaRt).")
}

bulk_mat2 <- bulk_mat[common_genes, , drop = FALSE]
sce2 <- sce[common_genes, ]

# ---------------------------
# 4) Helper: pick genes that distinguish cell types *within a group*
#    Improved with lfc_min threshold for better specificity
# ---------------------------
pick_within_group_markers <- function(sce_obj, types, clusters_col="segment_use",
                                      top_n=120, min_mean=0.05, lfc_min=1.0) {
  X <- assay(sce_obj, "counts")
  cl <- as.character(colData(sce_obj)[[clusters_col]])

  keep <- cl %in% types
  if (sum(keep) == 0) return(character(0))
  
  Xg <- X[, keep, drop=FALSE]
  clg <- factor(cl[keep])

  # mean expression per type
  means <- sapply(levels(clg), function(t) Matrix::rowMeans(Xg[, clg == t, drop=FALSE]))
  rownames(means) <- rownames(Xg)

  markers <- character(0)
  for (t in colnames(means)) {
    mu_t <- means[, t]
    mu_o <- rowMeans(means[, setdiff(colnames(means), t), drop=FALSE])
    lfc  <- log2((mu_t + 1) / (mu_o + 1))

    # Stricter filtering: require minimum LFC (2x = lfc 1.0)
    ok <- is.finite(lfc) & (mu_t > min_mean) & (lfc > lfc_min)
    if (sum(ok) == 0) next

    top <- names(sort(lfc[ok], decreasing=TRUE))[seq_len(min(top_n, sum(ok)))]
    markers <- unique(c(markers, top))
  }
  markers
}

message("\n--- Building within-group marker sets ---")
group.markers <- list(
  Glomerular = pick_within_group_markers(sce2, clusters.type$Glomerular, top_n=120, min_mean=0.05, lfc_min=1.0),
  Tubule     = pick_within_group_markers(sce2, clusters.type$Tubule,     top_n=120, min_mean=0.05, lfc_min=1.0),
  Immune     = NULL,
  Stroma     = NULL
)

message("Glomerular markers: ", length(group.markers$Glomerular))
message("Tubule markers: ", length(group.markers$Tubule))

# MuSiC expects genes x samples numeric matrix
stopifnot(is.matrix(bulk_mat2))
stopifnot(is.numeric(bulk_mat2))
stopifnot(!is.null(rownames(bulk_mat2)))
stopifnot(!is.null(colnames(bulk_mat2)))

# ---------------------------
# 5) Run tree-guided MuSiC (music_prop.cluster)
# ---------------------------
message("\nRunning tree-guided MuSiC deconvolution (music_prop.cluster)...")

stopifnot("segment_use" %in% colnames(colData(sce2)))
stopifnot("clusterType" %in% colnames(colData(sce2)))
stopifnot("donor_id" %in% colnames(colData(sce2)))

# Ensure character columns for MuSiC
colData(sce2)$segment_use <- as.character(colData(sce2)$segment_use)
colData(sce2)$clusterType <- as.character(colData(sce2)$clusterType)
colData(sce2)$donor_id <- as.character(colData(sce2)$donor_id)

res_cluster <- MuSiC::music_prop.cluster(
  bulk.mtx      = bulk_mat2,        # genes x samples
  sc.sce        = sce2,
  group.markers = group.markers,
  groups        = "clusterType",
  clusters      = "segment_use",
  samples       = "donor_id",       # must exist in colData(sce2)
  clusters.type = clusters.type,
  normalize     = FALSE,            # IMPORTANT (see MuSiC docs)
  centered      = FALSE,
  verbose       = TRUE
)

# Extract cell-type proportions
prop_seg <- res_cluster$Est.prop.weighted.cluster
write.csv(prop_seg, file.path(outdir, "music_cluster_proportions.csv"), quote = FALSE)
message("Wrote tree-guided cluster proportions to: ", file.path(outdir, "music_cluster_proportions.csv"))

# Compute coarse group proportions by summing segments (MuSiC doesn't return this directly)
prop_group <- data.frame(
  Glomerular = prop_seg[, "Podocyte"] + prop_seg[, "Endothelial"] + prop_seg[, "Mesangial"],
  Tubule = prop_seg[, "PT"] + prop_seg[, "TAL_LOH"] + prop_seg[, "DCT"] + 
           prop_seg[, "CD"] + prop_seg[, "Other"],
  Immune = prop_seg[, "Immune"],
  Stroma = prop_seg[, "Fibroblast"]
)
rownames(prop_group) <- rownames(prop_seg)
write.csv(prop_group, file.path(outdir, "music_group_proportions.csv"), quote = FALSE)
message("Wrote coarse group proportions to: ", file.path(outdir, "music_group_proportions.csv"))

# Quick summary of all segments
message("\n--- Segment Proportion Summary ---")
print(round(colMeans(prop_seg), 4))

message("\n--- Coarse Group Proportion Summary ---")
print(round(colMeans(prop_group), 4))

# ---------------------------
# 5b) Generalized Validation: Data-Driven Markers for ALL Cell Types
# ---------------------------
message("\n", paste(rep("=", 60), collapse = ""))
message("GENERALIZED VALIDATION: All Segments")
message(paste(rep("=", 60), collapse = ""))

# Function to get top markers per segment from the reference (data-driven)
get_top_markers_per_segment <- function(sce_obj, clusters_col = "segment_use", 
                                         top_n = 20, min_mean = 0.05, min_lfc = 1.0) {
  X <- assay(sce_obj, "counts")
  cl <- as.character(colData(sce_obj)[[clusters_col]])
  
  segs <- unique(cl)
  # Mean expression per segment
  means <- sapply(segs, function(s) Matrix::rowMeans(X[, cl == s, drop = FALSE]))
  colnames(means) <- segs
  
  # For each segment, find genes with highest fold-change vs all others
markers <- lapply(segs, function(s) {
    mu_s <- means[, s]
    mu_other <- rowMeans(means[, setdiff(segs, s), drop = FALSE])
    lfc <- log2((mu_s + 1) / (mu_other + 1))
    
    ok <- is.finite(lfc) & (mu_s > min_mean) & (lfc > min_lfc)
    if (sum(ok) == 0) return(character(0))
    
    names(sort(lfc[ok], decreasing = TRUE))[seq_len(min(top_n, sum(ok)))]
  })
  names(markers) <- segs
  markers
}

# Get top 20 markers for each segment from reference
message("\nDeriving top markers per segment from reference...")
segment_markers <- get_top_markers_per_segment(sce2, top_n = 20, min_lfc = 1.0)

message("\nMarkers found per segment:")
for (seg in names(segment_markers)) {
  message(sprintf("  %s: %d markers", seg, length(segment_markers[[seg]])))
}

# ---------------------------
# FIX APPLIED 2026-01-11: CPM Normalization for Validation
# ---------------------------
# Problem: prop_seg is library-size independent (fractions sum to 1), but raw 
#          bulk counts scale with sequencing depth. This caused spurious low/negative
#          correlations in validation (e.g., DCT showed ρ=0.03 with raw, ρ=0.33 with CPM).
# Solution: Normalize bulk to CPM before computing marker scores.
# ---------------------------

# CPM normalize bulk matrix for validation
lib_sizes <- colSums(bulk_mat2)
bulk_cpm <- sweep(bulk_mat2, 2, lib_sizes, "/") * 1e6
message("CPM normalized bulk matrix for validation (library sizes: ", 
        round(min(lib_sizes)/1e6, 1), "M - ", round(max(lib_sizes)/1e6, 1), "M)")

# Function to compute bulk expression scores using each segment's markers (CPM-normalized)
compute_bulk_scores <- function(bulk_cpm_mat, markers_list) {
  scores <- sapply(names(markers_list), function(seg) {
    genes <- intersect(markers_list[[seg]], rownames(bulk_cpm_mat))
    if (length(genes) == 0) return(rep(NA, ncol(bulk_cpm_mat)))
    colMeans(log1p(bulk_cpm_mat[genes, , drop = FALSE]))
  })
  as.data.frame(scores)
}

bulk_scores <- compute_bulk_scores(bulk_cpm, segment_markers)


# Validate: correlate estimated proportions with bulk marker scores
message("\n--- Validation: Proportion vs Bulk Marker Score Correlations ---")

validation_results <- data.frame(
  segment = colnames(prop_seg),
  n_markers = sapply(colnames(prop_seg), function(s) {
    if (s %in% names(segment_markers)) {
      length(intersect(segment_markers[[s]], rownames(bulk_mat2)))
    } else NA
  }),
  spearman_rho = NA,
  pvalue = NA,
  interpretation = NA
)

for (i in seq_len(nrow(validation_results))) {
  seg <- validation_results$segment[i]
  if (!(seg %in% colnames(bulk_scores))) next
  if (all(is.na(bulk_scores[[seg]]))) next
  
  ct <- cor.test(prop_seg[, seg], bulk_scores[[seg]], method = "spearman")
  validation_results$spearman_rho[i] <- round(ct$estimate, 3)
  validation_results$pvalue[i] <- signif(ct$p.value, 3)
  
  # Interpretation
  rho <- ct$estimate
  validation_results$interpretation[i] <- if (rho > 0.5) {
    "Strong"
  } else if (rho > 0.3) {
    "Moderate"
  } else if (rho > 0.1) {
    "Weak"
  } else if (rho > -0.1) {
    "None"
  } else {
    "Negative!"
  }
}

print(validation_results)

# Save validation results
write.csv(validation_results, file.path(outdir, "validation_segment_correlations.csv"), 
          row.names = FALSE)
message("\nWrote validation results to: ", file.path(outdir, "validation_segment_correlations.csv"))

# Summary statistics
good_segments <- validation_results$segment[!is.na(validation_results$spearman_rho) & 
                                              validation_results$spearman_rho > 0.3]
poor_segments <- validation_results$segment[!is.na(validation_results$spearman_rho) & 
                                              validation_results$spearman_rho < 0.1]

message("\n--- Validation Summary ---")
message("Segments with good correlation (rho > 0.3): ", 
        if (length(good_segments) > 0) paste(good_segments, collapse = ", ") else "None")
message("Segments with poor correlation (rho < 0.1): ", 
        if (length(poor_segments) > 0) paste(poor_segments, collapse = ", ") else "None")

# ---------------------------
# 5c) Coarse Group Validation
# ---------------------------
message("\n--- Coarse Group Split Sanity Check ---")

for (grp in names(clusters.type)) {
  segs <- intersect(clusters.type[[grp]], colnames(prop_seg))
  if (length(segs) == 0) next
  
  grp_total <- rowSums(prop_seg[, segs, drop = FALSE])
  message(sprintf("\n%s (segments: %s):", grp, paste(segs, collapse = ", ")))
  message(sprintf("  Mean total: %.3f, Median: %.3f, Range: [%.4f, %.3f]",
                  mean(grp_total), median(grp_total), min(grp_total), max(grp_total)))
  
  # Within-group breakdown
  if (length(segs) > 1) {
    for (s in segs) {
      frac <- prop_seg[, s] / (grp_total + 1e-9)
      message(sprintf("    %s fraction within %s: median=%.2f", s, grp, median(frac)))
    }
  }
}


























# ---------------------------
# 6) Use segment-direct proportions (already at segment level)
# ---------------------------
# Since we ran MuSiC on segments directly, prop_seg is already at segment granularity
# No need to aggregate cell types - the output is already segment-level

segment_prop <- as.data.frame(prop_seg)

message("\n--- Segment Proportions Summary ---")
print(colnames(segment_prop))
print(colMeans(segment_prop, na.rm = TRUE))

# ---------------------------
# 7) CLR transform for downstream regression (compositional handling)
# ---------------------------
# Your pipeline explicitly suggests CLR (or dropping one category).

clr_transform <- function(P, pseudo = 1e-6) {
  P2 <- pmax(P, pseudo)
  logP <- log(P2)
  gm <- rowMeans(logP)
  sweep(logP, 1, gm, "-")
}

segment_clr <- clr_transform(as.matrix(segment_prop))
write.csv(segment_clr, file.path(outdir, "music_segment_direct_proportions_CLR.csv"), quote = FALSE)
message("Wrote CLR-transformed segment proportions to: ", file.path(outdir, "music_segment_direct_proportions_CLR.csv"))

# ---------------------------
# 8) Quick QC plots
# ---------------------------
df_long <- data.frame(
  sample = rownames(segment_prop),
  segment_prop,
  check.names = FALSE
)
df_melt <- data.table::melt(
  data.table::as.data.table(df_long),
  id.vars = "sample",
  variable.name = "segment",
  value.name = "prop"
)
df_melt <- as.data.frame(df_melt)


p <- ggplot(df_melt, aes(x = segment, y = prop)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MuSiC segment proportion distributions", y = "Estimated proportion")
ggsave(file.path(outdir, "qc_segment_prop_boxplot.png"), p, width = 9, height = 4)
message("Saved QC plot to: ", file.path(outdir, "qc_segment_prop_boxplot.png"))

message("Done. Wrote outputs to: ", normalizePath(outdir))
