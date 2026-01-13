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
bulk_meta_csv <- "data/processed/aligned_outputs/metadata_aligned.tsv"

# Use the H5AD you wrote (best), NOT the obs/var CSV
# e.g. tms_kidney_female/tms_kidney_female_ALLDATASETS_innerGenes.h5ad
sc_h5ad_path <- "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"


outdir <- "data/processed/deconvolution/test16"
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
  if (length(ids) == 0) {
    return("unknown")
  }
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
obs <- data.frame(fread(obs_path, header = TRUE), row.names = 1, check.names = FALSE)
var <- data.frame(fread(var_path, header = TRUE), row.names = 1, check.names = FALSE)

# Handle potential transpose issues with mtx (Scanpy writes shape n_obs x n_vars, readMM reads as such)
# But SCE expects genes x cells (n_vars x n_obs). Scanpy usually writes cells x genes?
# Creating scipy.io.mmwrite(..., raw_adata.X) usually writes standard MM.
# If raw_adata.X is (cells, genes), then readMM returns (rows=cells, cols=genes).
# SCE needs (rows=genes, cols=cells). So we likely need to Transpose.
# Let's check dimensions carefully.
message("MTX dims: ", paste(dim(raw_mat), collapse = " x "))
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
  if (is.na(x)) {
    return(NA_character_)
  }
  x2 <- tolower(trimws(x))

  if (grepl("podocyte", x2)) {
    return("Podocyte")
  }
  if (grepl("mesangial", x2)) {
    return("Mesangial")
  }
  if (grepl("pecam|endothelial|capillary|artery", x2)) {
    return("Endothelial")
  }
  if (grepl("fibroblast|stroma|stromal", x2)) {
    return("Fibroblast")
  }
  if (grepl("^cd45|macrophage|\\bt\\s*cell\\b|\\bb\\s*cell\\b|plasma\\s*cell|nk\\s*cell|lymph|leukocyte", x2)) {
    return("Immune")
  }

  # Tubules (allow words in-between)
  if (grepl("proximal.*tubule|proximal\\s+tube|\\bpt\\b", x2)) {
    return("PT")
  }
  if (grepl("distal.*convoluted|distal.*tubule|\\bdct\\b", x2)) {
    return("DCT")
  }
  if (grepl("thick.*ascending|loop\\s+of\\s+henle|henle|\\btal\\b", x2)) {
    return("TAL_LOH")
  }
  if (grepl("collecting\\s+duct|principal\\s+cell|intercalated", x2)) {
    return("CD")
  }

  if (grepl("brush\\s*cell|tuft\\s*cell", x2)) {
    return("Other")
  }
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
# 2) REDEFINED Coarse groups (Plan B: Split Proximal vs Distal)
# ----------------------------
seg <- as.character(colData(sce)$segment_use)

# Define groupings: Separate PT from the Distal Nephron
colData(sce)$clusterType <- ifelse(seg %in% c("Podocyte", "Endothelial", "Mesangial"), "Glomerular",
  ifelse(seg %in% c("Immune"), "Immune",
    ifelse(seg %in% c("Fibroblast"), "Stroma",
      ifelse(seg %in% c("PT"), "Proximal", "Distal") # <--- THE CRITICAL SPLIT
    )
  )
)
colData(sce)$clusterType <- factor(colData(sce)$clusterType)

# Define the tree structure for MuSiC
clusters.type <- list(
  Glomerular = c("Podocyte", "Endothelial", "Mesangial"),
  Proximal   = c("PT"), # PT stands alone
  Distal     = c("TAL_LOH", "DCT", "CD", "Other"), # Distal forms a coalition
  Immune     = c("Immune"),
  Stroma     = c("Fibroblast")
)

message("\n--- New Coarse Group Distribution (Proximal vs Distal Split) ---")
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
# EXPERIMENTAL FLAGS (set TRUE to enable diagnostics)
# ---------------------------
EXP_REMOVE_OTHER <- TRUE # Experiment A: remove "Other" from Tubule group
EXP_PER_DONOR_PSEUDOBULK <- FALSE # Experiment B: sample from single donor per mixture
EXP_USE_CANONICAL_MARKERS <- TRUE # Experiment C: use canonical tubule markers

# ---------------------------
# 4) Helper: pick genes that distinguish cell types *within a group*
#    FIXED: use logCPM for marker selection (raw counts are depth-confounded)
# ---------------------------
pick_within_group_markers <- function(sce_obj, types, clusters_col = "segment_use",
                                      top_n = 120, min_mean = 0.5, lfc_min = 1.0) {
  X <- assay(sce_obj, "counts")
  cl <- as.character(colData(sce_obj)[[clusters_col]])

  keep <- cl %in% types
  if (sum(keep) == 0) {
    return(character(0))
  }

  Xg <- X[, keep, drop = FALSE]
  clg <- factor(cl[keep])

  # FIX: compute logCPM for marker selection (removes depth confounding)
  lib <- Matrix::colSums(Xg)
  Xg_cpm <- t(t(Xg) / (lib + 1e-12)) * 1e6
  Xg_log <- log1p(Xg_cpm) # logCPM for normalized marker derivation

  # mean logCPM expression per type
  means <- sapply(levels(clg), function(t) Matrix::rowMeans(Xg_log[, clg == t, drop = FALSE]))
  rownames(means) <- rownames(Xg)

  markers <- character(0)
  for (t in colnames(means)) {
    mu_t <- means[, t]
    mu_o <- rowMeans(means[, setdiff(colnames(means), t), drop = FALSE])
    # LFC on logCPM scale (no pseudocount needed since log1p already applied)
    lfc <- mu_t - mu_o

    # min_mean now on logCPM scale (0.5 ≈ 0.65 CPM minimum)
    ok <- is.finite(lfc) & (mu_t > min_mean) & (lfc > lfc_min)
    if (sum(ok) == 0) next

    top <- names(sort(lfc[ok], decreasing = TRUE))[seq_len(min(top_n, sum(ok)))]
    markers <- unique(c(markers, top))
  }
  markers
}

# ---------------------------
# Helper: pick markers for a SINGLE target type vs other types
# Balanced marker selection per subtype
# ---------------------------
pick_type_markers <- function(sce_obj, target, types,
                              clusters_col = "segment_use",
                              top_n = 80, min_mean = 0.1, lfc_min = 0.5) {
  X <- counts(sce_obj)
  cl <- as.character(colData(sce_obj)[[clusters_col]])
  keep <- cl %in% types
  Xg <- X[, keep, drop = FALSE]
  clg <- factor(cl[keep])

  lib <- Matrix::colSums(Xg)
  Xcpm <- t(t(Xg) / (lib + 1e-12)) * 1e6
  Xlog <- log1p(Xcpm)

  mu <- sapply(levels(clg), function(tt) Matrix::rowMeans(Xlog[, clg == tt, drop = FALSE]))
  rownames(mu) <- rownames(Xg)

  mu_t <- mu[, target]
  mu_o <- rowMeans(mu[, setdiff(colnames(mu), target), drop = FALSE])
  lfc <- mu_t - mu_o

  ok <- is.finite(lfc) & (mu_t > min_mean) & (lfc > lfc_min)
  if (!any(ok)) {
    return(character(0))
  }

  head(names(sort(lfc[ok], decreasing = TRUE)), top_n)
}

# ---------------------------
# Experiment A: Remove "Other" from analysis
# ---------------------------
if (EXP_REMOVE_OTHER) {
  message("\n*** EXPERIMENT A: Removing 'Other' segment from analysis ***")
  keep_seg <- colData(sce2)$segment_use != "Other"
  sce2 <- sce2[, keep_seg]
  colData(sce2)$segment_use <- droplevels(factor(colData(sce2)$segment_use))

  # PLAN B: Rebuild clusterType WITH Proximal/Distal split (not single Tubule bucket)
  seg <- as.character(colData(sce2)$segment_use)
  colData(sce2)$clusterType <- factor(ifelse(
    seg %in% c("Podocyte", "Endothelial", "Mesangial"), "Glomerular",
    ifelse(seg == "Immune", "Immune",
      ifelse(seg == "Fibroblast", "Stroma",
        ifelse(seg == "PT", "Proximal", "Distal") # <--- PLAN B: PT isolated, TAL/DCT/CD -> Distal
      )
    )
  ))

  # PLAN B: clusters.type with Proximal/Distal split (no single Tubule bucket)
  clusters.type <- list(
    Glomerular = c("Podocyte", "Endothelial", "Mesangial"),
    Proximal   = c("PT"), # PT stands alone
    Distal     = c("TAL_LOH", "DCT", "CD"), # Distal coalition (Other removed)
    Immune     = c("Immune"),
    Stroma     = c("Fibroblast")
  )

  message("Segment distribution after removing 'Other':")
  print(table(colData(sce2)$segment_use))
  message("\n--- PLAN B: Proximal vs Distal Split Active ---")
  print(table(colData(sce2)$clusterType))
}

message("\n--- Building within-group marker sets (logCPM-based) ---")
# ---------------------------
# OPTIMIZED MARKER SELECTION (For Split Groups)
# ---------------------------
message("\n--- Building Optimized Marker List (Canonical Focus + Group Split) ---")
tub_types <- intersect(c("PT", "TAL_LOH", "DCT", "CD"), unique(colData(sce2)$segment_use))

# 1. PT Markers (For the "Proximal" Group)
mk_pt <- pick_type_markers(sce2, "PT", tub_types, top_n = 50)

# 2. Distal Markers (We build specific lists for the sub-types)
# TAL: Strict Canonical
mk_tal_data <- pick_type_markers(sce2, "TAL_LOH", tub_types, top_n = 15)
canon_tal <- c("Umod", "Slc12a1", "Cldn10", "Cldn16", "Kcnj1", "Bsnd")
canon_tal_res <- resolve_markers_to_matrix(canon_tal, rownames(sce2))
mk_tal_final <- unique(c(canon_tal_res, mk_tal_data))
message("Resolved ", length(canon_tal_res), " canonical TAL markers.")

# DCT: Strict Canonical
mk_dct_data <- pick_type_markers(sce2, "DCT", tub_types, top_n = 15)
canon_dct <- c("Slc12a3", "Pvalb", "Calb1", "Trpv5", "Kl", "Wnk4")
canon_dct_res <- resolve_markers_to_matrix(canon_dct, rownames(sce2))
mk_dct_final <- unique(c(canon_dct_res, mk_dct_data))
message("Resolved ", length(canon_dct_res), " canonical DCT markers.")

# CD: Strict Canonical
mk_cd_data <- pick_type_markers(sce2, "CD", tub_types, top_n = 15)
canon_cd <- c("Aqp2", "Aqp3", "Fxyd4", "Hsd11b2", "Scnn1g", "Krt8", "Krt18")
canon_cd_res <- resolve_markers_to_matrix(canon_cd, rownames(sce2))
mk_cd_final <- unique(c(canon_cd_res, mk_cd_data))
message("Resolved ", length(canon_cd_res), " canonical CD markers.")

# 3. PLAN B: Always use Proximal/Distal split structure
# This ensures PT variance is isolated from Distal subtypes (TAL/DCT/CD)
group.markers <- list(
  Glomerular = pick_within_group_markers(sce2, clusters.type$Glomerular, top_n = 50),
  Proximal   = mk_pt, # PT markers for Proximal group
  Distal     = unique(c(mk_tal_final, mk_dct_final, mk_cd_final)), # Combined Distal markers
  Immune     = NULL,
  Stroma     = NULL
)
message("Final Marker Counts (PLAN B - Proximal/Distal Split):")
message("  Proximal Group: ", length(group.markers$Proximal))
message("  Distal Group:   ", length(group.markers$Distal))

# ---------------------------
# Experiment C: Use canonical tubule markers (override data-driven)
# PLAN B: Apply separately to Proximal and Distal groups
# ---------------------------
if (EXP_USE_CANONICAL_MARKERS) {
  message("\n*** EXPERIMENT C: Using canonical tubule markers (PLAN B structure) ***")

  # Canonical markers for PROXIMAL (PT)
  canonical_proximal <- c("Slc34a1", "Lrp2", "Aqp1", "Slc22a6", "Slc5a2")
  resolved_proximal <- resolve_markers_to_matrix(canonical_proximal, rownames(sce2))
  message("Resolved ", length(resolved_proximal), " of ", length(canonical_proximal), " canonical Proximal (PT) markers")
  if (length(resolved_proximal) > 0) {
    group.markers$Proximal <- unique(c(resolved_proximal, mk_pt))
    message("Proximal canonical markers: ", paste(resolved_proximal, collapse = ", "))
  }

  # Canonical markers for DISTAL (TAL, DCT, CD)
  canonical_distal <- c(
    # TAL markers
    "Umod", "Slc12a1",
    # DCT markers
    "Slc12a3", "Pvalb", "Calb1",
    # CD markers
    "Aqp2", "Krt8", "Krt18", "Atp6v1b1"
  )
  resolved_distal <- resolve_markers_to_matrix(canonical_distal, rownames(sce2))
  message("Resolved ", length(resolved_distal), " of ", length(canonical_distal), " canonical Distal markers")
  if (length(resolved_distal) > 0) {
    group.markers$Distal <- unique(c(resolved_distal, mk_tal_final, mk_dct_final, mk_cd_final))
    message("Distal canonical markers: ", paste(resolved_distal, collapse = ", "))
  }
}

message("\n--- PLAN B Marker Summary ---")
message("Glomerular markers: ", length(group.markers$Glomerular))
message("Proximal markers: ", length(group.markers$Proximal))
message("Distal markers: ", length(group.markers$Distal))

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

# ---------------------------
# CRITICAL: donor_id diagnostics
# MuSiC uses cross-donor variance for cell-size (theta) estimation.
# If donor_id is NA, constant, or unique per cell, resolution fails.
# ---------------------------
message("\n", paste(rep("=", 60), collapse = ""))
message("DONOR_ID DIAGNOSTICS (MuSiC requires valid multi-donor data)")
message(paste(rep("=", 60), collapse = ""))

did <- colData(sce2)$donor_id
cat("donor_id NAs:", sum(is.na(did)), "of", length(did), "\n")

# Remove NAs for distribution analysis
did_clean <- did[!is.na(did)]
tab <- sort(table(did_clean), decreasing = TRUE)
cat("\nn donors:", length(tab), "\n")
cat("median cells/donor:", median(tab), "\n")
cat("min cells/donor:", min(tab), "\n")
cat("max cells/donor:", max(tab), "\n")

message("\nTop 20 donors by cell count:")
print(head(tab, 20))

# Check if each donor has cells in each segment
message("\nDonor x Segment coverage (does each donor have each segment?):")
cross_tab <- table(colData(sce2)$donor_id, colData(sce2)$segment_use)
coverage_mat <- cross_tab > 0
print(coverage_mat)

# Summary: how many segments per donor?
segs_per_donor <- rowSums(coverage_mat)
cat("\nSegments covered per donor:\n")
print(summary(segs_per_donor))

# Flag potential issues
n_donors <- length(tab)
if (n_donors == 1) {
  message("\n*** WARNING: Only 1 donor! MuSiC cell-size estimation will be unstable. ***")
  message("*** Consider: colData(sce2)$donor_id <- 'D1' to force single-donor mode ***")
} else if (n_donors < 5) {
  message("\n*** WARNING: Only ", n_donors, " donors. MuSiC may have limited variance estimation. ***")
} else if (n_donors > 0.5 * length(did_clean)) {
  message("\n*** WARNING: donor_id looks like barcodes (", n_donors, " 'donors' for ", length(did_clean), " cells) ***")
  message("*** This breaks MuSiC's cross-subject variance model! ***")
} else {
  message("\n✓ donor_id looks reasonable for MuSiC")
}

# Option to force single donor for testing
# Uncomment the line below to bypass donor issues:
# colData(sce2)$donor_id <- "D1"

message(paste(rep("=", 60), collapse = ""))

# Must include PT
print(table(colData(sce2)$segment_use))

# Must match clusters.type names
print(names(table(colData(sce2)$segment_use)))
print(clusters.type$Tubule)

# ============================================================
# In-silico pseudo-bulk test for MuSiC
# Requires: sce2 with counts, colData(sce2)$segment_use and $clusterType and $donor_id
# Also requires: clusters.type and group.markers as in your script
# ============================================================

set.seed(1)

# ---------
# Helpers
# ---------
normalize_rows_to_one <- function(M, eps = 1e-12) {
  rs <- rowSums(M)
  rs[rs < eps] <- 1
  M / rs
}

rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))

# Dirichlet sampler (no extra packages)
rdirichlet_base <- function(n, alpha) {
  k <- length(alpha)
  X <- matrix(0, nrow = n, ncol = k)
  for (j in seq_len(k)) X[, j] <- rgamma(n, shape = alpha[j], rate = 1)
  X <- X / rowSums(X)
  X
}

# ---------
# IMPROVED: Pseudobulk with exact proportions and all fixes
# ---------
make_pseudobulk_exact <- function(
  sce_obj,
  segment_col = "segment_use",
  donor_col = "donor_id",
  n_mixtures = 40,
  cells_per_mixture = 1500,
  depth_scale = 1.0,
  target_depth = NULL, # Fixed depth (e.g., 20e6), overrides depth_scale
  use_multinomial_depth = TRUE, # Use multinomial downsampling instead of Poisson
  # Two-level design: independent ranges for tubule segments
  dct_range = c(0.00, 0.25),
  pt_range = c(0.40, 0.75),
  tal_range = c(0.05, 0.25),
  cd_range = c(0.00, 0.15),
  other_tubule_range = c(0.00, 0.10),
  # Coarse group allocation
  tubule_total_range = c(0.50, 0.85),
  glom_total_range = c(0.05, 0.20),
  immune_frac_range = c(0.01, 0.10),
  stroma_frac_range = c(0.02, 0.10),
  seed = 1
) {
  set.seed(seed)
  stopifnot("counts" %in% assayNames(sce_obj))

  seg <- factor(as.character(colData(sce_obj)[[segment_col]]))
  donor <- as.character(colData(sce_obj)[[donor_col]])
  seg_levels <- levels(seg)
  donors <- unique(donor[!is.na(donor)])

  genes <- rownames(sce_obj)
  n_genes <- length(genes)

  # Define coarse groups based on segment names
  glom_segs <- intersect(c("Podocyte", "Endothelial", "Mesangial"), seg_levels)
  tubule_segs <- intersect(c("PT", "TAL_LOH", "DCT", "CD", "Other"), seg_levels)
  immune_segs <- intersect(c("Immune"), seg_levels)
  stroma_segs <- intersect(c("Fibroblast"), seg_levels)

  # Pre-compute index by donor x segment
  cell_idx <- seq_len(ncol(sce_obj))
  idx_by_donor_seg <- list()
  for (d in donors) {
    idx_by_donor_seg[[d]] <- split(cell_idx[donor == d], seg[donor == d])
  }

  # Check how many segments each donor covers
  seg_coverage <- sapply(donors, function(d) sum(sapply(idx_by_donor_seg[[d]], length) > 0))
  good_donors <- donors[seg_coverage >= 4] # Need at least 4 segments for reasonable mixtures
  if (length(good_donors) == 0) good_donors <- donors

  message("Pseudobulk exact: ", length(good_donors), " donors available with >=4 segments")

  # Output matrices
  pb_mat <- matrix(0, nrow = n_genes, ncol = n_mixtures, dimnames = list(genes, paste0("pb", seq_len(n_mixtures))))
  P_target <- matrix(0, nrow = n_mixtures, ncol = length(seg_levels), dimnames = list(paste0("pb", seq_len(n_mixtures)), seg_levels))
  P_cell_true <- matrix(0, nrow = n_mixtures, ncol = length(seg_levels), dimnames = list(paste0("pb", seq_len(n_mixtures)), seg_levels))
  P_rna_true <- matrix(0, nrow = n_mixtures, ncol = length(seg_levels), dimnames = list(paste0("pb", seq_len(n_mixtures)), seg_levels))
  used_counts <- matrix(0, nrow = n_mixtures, ncol = length(seg_levels), dimnames = list(paste0("pb", seq_len(n_mixtures)), seg_levels))
  donor_used <- character(n_mixtures)

  for (i in seq_len(n_mixtures)) {
    # ==========================================
    # FIX 5: Per-donor sampling
    # Pick one donor for this mixture
    # ==========================================
    d <- sample(good_donors, 1)
    donor_used[i] <- d
    idx_this_donor <- idx_by_donor_seg[[d]]

    # ==========================================
    # FIX 3: Independent DCT variation via two-level design
    # Sample tubule proportions with independent ranges
    # ==========================================

    # Step 1: Sample coarse group totals
    tubule_total <- runif(1, tubule_total_range[1], tubule_total_range[2])
    glom_total <- runif(1, glom_total_range[1], glom_total_range[2])
    immune_total <- runif(1, immune_frac_range[1], immune_frac_range[2])
    stroma_total <- runif(1, stroma_frac_range[1], stroma_frac_range[2])

    # Normalize coarse groups to sum to 1
    coarse_raw <- c(Tubule = tubule_total, Glomerular = glom_total, Immune = immune_total, Stroma = stroma_total)
    coarse_norm <- coarse_raw / sum(coarse_raw)

    # Step 2: Within-tubule, sample each segment independently with custom ranges
    tubule_within <- c(
      DCT = runif(1, dct_range[1], dct_range[2]),
      PT = runif(1, pt_range[1], pt_range[2]),
      TAL_LOH = runif(1, tal_range[1], tal_range[2]),
      CD = runif(1, cd_range[1], cd_range[2]),
      Other = runif(1, other_tubule_range[1], other_tubule_range[2])
    )

    # Keep only segments that exist in this dataset
    tubule_within <- tubule_within[names(tubule_within) %in% tubule_segs]
    if (length(tubule_within) > 0) {
      tubule_within <- tubule_within / sum(tubule_within) # normalize within tubule
    }

    # Step 3: Within-glomerular, equal allocation or Dirichlet
    if (length(glom_segs) > 0) {
      glom_within <- rdirichlet_base(1, rep(2, length(glom_segs)))[1, ]
      names(glom_within) <- glom_segs
    } else {
      glom_within <- numeric(0)
    }

    # Step 4: Combine into final segment target proportions
    target_props <- setNames(rep(0, length(seg_levels)), seg_levels)

    for (s in names(tubule_within)) {
      if (s %in% seg_levels) target_props[s] <- coarse_norm["Tubule"] * tubule_within[s]
    }
    for (s in names(glom_within)) {
      if (s %in% seg_levels) target_props[s] <- coarse_norm["Glomerular"] * glom_within[s]
    }
    for (s in immune_segs) {
      if (s %in% seg_levels) target_props[s] <- coarse_norm["Immune"] / length(immune_segs)
    }
    for (s in stroma_segs) {
      if (s %in% seg_levels) target_props[s] <- coarse_norm["Stroma"] / length(stroma_segs)
    }

    # Normalize target proportions to sum to 1
    if (sum(target_props) > 0) {
      target_props <- target_props / sum(target_props)
    }

    P_target[i, ] <- target_props

    # ==========================================
    # FIX 2: Multinomial for exact cell counts
    # P_target is the soft target, rmultinom gives exact counts
    # ==========================================
    n_per_seg <- as.integer(rmultinom(1, size = cells_per_mixture, prob = target_props))
    names(n_per_seg) <- seg_levels
    used_counts[i, ] <- n_per_seg

    # ==========================================
    # Sample cells and track UMIs per segment
    # ==========================================
    umi_by_seg <- setNames(rep(0, length(seg_levels)), seg_levels)
    x_total <- numeric(n_genes)
    names(x_total) <- genes

    for (s in seg_levels) {
      if (n_per_seg[s] <= 0) next
      pool <- idx_this_donor[[s]]
      if (is.null(pool) || length(pool) == 0) {
        # This donor doesn't have this segment, skip
        next
      }

      chosen_s <- sample(pool, size = n_per_seg[s], replace = TRUE)
      xs <- Matrix::rowSums(counts(sce_obj)[, chosen_s, drop = FALSE])
      xs <- as.numeric(xs)
      names(xs) <- genes

      umi_by_seg[s] <- sum(xs)
      x_total <- x_total + xs
    }

    # ==========================================
    # FIX 1: Compute P_rna_true from actual UMI contributions
    # ==========================================
    total_umi <- sum(umi_by_seg)
    if (total_umi > 0) {
      P_rna_true[i, ] <- umi_by_seg / total_umi
    }

    # P_cell_true from actual cell counts (adjusted for any missing segments)
    actual_cells <- used_counts[i, ]
    # Adjust for cells we couldn't sample (donor didn't have that segment)
    for (s in seg_levels) {
      if (umi_by_seg[s] == 0 && used_counts[i, s] > 0) {
        actual_cells[s] <- 0 # Couldn't actually get these cells
      }
    }
    if (sum(actual_cells) > 0) {
      P_cell_true[i, ] <- actual_cells / sum(actual_cells)
    }

    # ==========================================
    # FIX 4: Multinomial downsampling for library size
    # Instead of per-gene Poisson noise, downsample to target depth
    # ==========================================
    if (use_multinomial_depth && sum(x_total) > 0) {
      actual_depth <- sum(x_total)
      if (!is.null(target_depth)) {
        final_depth <- min(target_depth, actual_depth) # Can't upsample
      } else {
        final_depth <- round(actual_depth * depth_scale)
      }

      # Multinomial downsampling
      prob_vec <- x_total / sum(x_total)
      x_total <- as.vector(rmultinom(1, size = final_depth, prob = prob_vec))
      names(x_total) <- genes
    }

    pb_mat[, i] <- x_total
  }

  message("Pseudobulk exact: generated ", n_mixtures, " mixtures from ", length(unique(donor_used)), " unique donors")
  message("  DCT range in P_target: [", round(min(P_target[, "DCT"]), 4), ", ", round(max(P_target[, "DCT"]), 4), "]")
  message("  PT range in P_target: [", round(min(P_target[, "PT"]), 4), ", ", round(max(P_target[, "PT"]), 4), "]")

  list(
    pb_mat = pb_mat,
    P_target = P_target,
    P_cell_true = P_cell_true,
    P_rna_true = P_rna_true,
    used_counts = used_counts,
    donor_used = donor_used
  )
}

# ---------
# Build pseudo-bulk mixtures
# ---------
make_pseudobulk <- function(
  sce_obj,
  segment_col = "segment_use",
  n_mixtures = 30,
  cells_per_mixture = 800,
  alpha = NULL,
  min_cells_per_segment = 10,
  use_poisson_depth = TRUE,
  depth_scale = 1.0
) {
  stopifnot("counts" %in% assayNames(sce_obj))
  seg <- factor(as.character(colData(sce_obj)[[segment_col]]))
  seg_levels <- levels(seg)

  if (is.null(alpha)) {
    alpha <- rep(0.7, length(seg_levels))
  }
  names(alpha) <- seg_levels

  idx_by_seg <- split(seq_len(ncol(sce_obj)), seg)

  # target proportions
  P_draw <- rdirichlet_base(n_mixtures, alpha = alpha)
  colnames(P_draw) <- seg_levels
  rownames(P_draw) <- paste0("pb", seq_len(n_mixtures))

  genes <- rownames(sce_obj)
  pb_mat <- matrix(0,
    nrow = length(genes), ncol = n_mixtures,
    dimnames = list(genes, rownames(P_draw))
  )

  # Two ground-truth matrices:
  # P_cell: fraction of sampled cells from each segment
  # P_rna:  fraction of total UMIs contributed by each segment
  P_cell <- matrix(0,
    nrow = n_mixtures, ncol = length(seg_levels),
    dimnames = list(rownames(P_draw), seg_levels)
  )
  P_rna <- matrix(0,
    nrow = n_mixtures, ncol = length(seg_levels),
    dimnames = list(rownames(P_draw), seg_levels)
  )

  for (i in seq_len(n_mixtures)) {
    n_per_seg <- pmax(round(P_draw[i, ] * cells_per_mixture), 0)

    # enforce minimum only for segments that are present
    n_per_seg[n_per_seg > 0 & n_per_seg < min_cells_per_segment] <- min_cells_per_segment

    # rescale to exactly cells_per_mixture
    if (sum(n_per_seg) > 0) {
      n_per_seg <- round(n_per_seg / sum(n_per_seg) * cells_per_mixture)
    }

    # one more adjust to hit exact total (optional but nice)
    diff <- cells_per_mixture - sum(n_per_seg)
    if (diff != 0) {
      j <- which.max(n_per_seg) # dump remainder into largest segment
      n_per_seg[j] <- n_per_seg[j] + diff
      n_per_seg[n_per_seg < 0] <- 0
    }

    # Track UMIs and cell counts per segment
    umi_by_seg <- setNames(rep(0, length(seg_levels)), seg_levels)
    n_by_seg <- setNames(rep(0, length(seg_levels)), seg_levels)

    x_total <- Matrix::Matrix(0, nrow = nrow(counts(sce_obj)), ncol = 1, sparse = TRUE)

    for (s in seg_levels) {
      if (n_per_seg[s] <= 0) next
      pool <- idx_by_seg[[s]]
      if (length(pool) == 0) next

      chosen_s <- sample(pool, size = n_per_seg[s], replace = TRUE)
      xs <- Matrix::rowSums(counts(sce_obj)[, chosen_s, drop = FALSE])
      umi_by_seg[s] <- sum(xs)
      n_by_seg[s] <- length(chosen_s)

      x_total <- x_total + xs
    }

    x <- as.numeric(x_total)
    names(x) <- genes

    if (use_poisson_depth) {
      lam <- x * depth_scale
      x <- rpois(length(lam), lambda = lam)
      names(x) <- genes
    }

    pb_mat[, i] <- x

    P_cell[i, ] <- n_by_seg / sum(n_by_seg)
    P_rna[i, ] <- umi_by_seg / sum(umi_by_seg)
  }

  list(pb_mat = pb_mat, P_cell = P_cell, P_rna = P_rna)
}

# ---------
# Experiment B: Per-donor pseudo-bulk (sample from single donor per mixture)
# ---------
make_pseudobulk_per_donor <- function(
  sce_obj,
  segment_col = "segment_use",
  donor_col = "donor_id",
  n_mixtures = 30,
  cells_per_mixture = 800,
  alpha = NULL,
  min_cells_per_segment = 5,
  use_poisson_depth = TRUE,
  depth_scale = 1.0
) {
  stopifnot("counts" %in% assayNames(sce_obj))
  seg <- factor(as.character(colData(sce_obj)[[segment_col]]))
  donor <- as.character(colData(sce_obj)[[donor_col]])
  seg_levels <- levels(seg)
  donors <- unique(donor[!is.na(donor)])

  if (length(donors) < 2) {
    message("WARNING: Only ", length(donors), " donor(s). Per-donor sampling not meaningful, falling back to pooled.")
    donors <- donors[1]
  }

  if (is.null(alpha)) {
    alpha <- rep(0.7, length(seg_levels))
  }
  names(alpha) <- seg_levels

  # Pre-compute index by donor x segment
  cell_idx <- seq_len(ncol(sce_obj))
  idx_by_donor_seg <- list()
  for (d in donors) {
    idx_by_donor_seg[[d]] <- split(cell_idx[donor == d], seg[donor == d])
  }

  # Check how many segments each donor has
  seg_coverage <- sapply(donors, function(d) sum(sapply(idx_by_donor_seg[[d]], length) > 0))
  good_donors <- donors[seg_coverage >= 3] # Need at least 3 segments
  if (length(good_donors) == 0) good_donors <- donors

  P_draw <- rdirichlet_base(n_mixtures, alpha = alpha)
  colnames(P_draw) <- seg_levels
  rownames(P_draw) <- paste0("pb", seq_len(n_mixtures))

  genes <- rownames(sce_obj)
  pb_mat <- matrix(0,
    nrow = length(genes), ncol = n_mixtures,
    dimnames = list(genes, rownames(P_draw))
  )

  P_cell <- matrix(0,
    nrow = n_mixtures, ncol = length(seg_levels),
    dimnames = list(rownames(P_draw), seg_levels)
  )
  P_rna <- matrix(0,
    nrow = n_mixtures, ncol = length(seg_levels),
    dimnames = list(rownames(P_draw), seg_levels)
  )
  donor_used <- character(n_mixtures)

  for (i in seq_len(n_mixtures)) {
    # Pick a random donor for this mixture
    d <- sample(good_donors, 1)
    donor_used[i] <- d
    idx_this <- idx_by_donor_seg[[d]]

    n_per_seg <- pmax(round(P_draw[i, ] * cells_per_mixture), 0)
    n_per_seg[n_per_seg > 0 & n_per_seg < min_cells_per_segment] <- min_cells_per_segment
    if (sum(n_per_seg) > 0) {
      n_per_seg <- round(n_per_seg / sum(n_per_seg) * cells_per_mixture)
    }
    diff <- cells_per_mixture - sum(n_per_seg)
    if (diff != 0) {
      j <- which.max(n_per_seg)
      n_per_seg[j] <- n_per_seg[j] + diff
      n_per_seg[n_per_seg < 0] <- 0
    }

    umi_by_seg <- setNames(rep(0, length(seg_levels)), seg_levels)
    n_by_seg <- setNames(rep(0, length(seg_levels)), seg_levels)
    x_total <- Matrix::Matrix(0, nrow = nrow(counts(sce_obj)), ncol = 1, sparse = TRUE)

    for (s in seg_levels) {
      if (n_per_seg[s] <= 0) next
      pool <- idx_this[[s]]
      if (is.null(pool) || length(pool) == 0) next

      chosen_s <- sample(pool, size = min(n_per_seg[s], length(pool)), replace = TRUE)
      xs <- Matrix::rowSums(counts(sce_obj)[, chosen_s, drop = FALSE])
      umi_by_seg[s] <- sum(xs)
      n_by_seg[s] <- length(chosen_s)
      x_total <- x_total + xs
    }

    x <- as.numeric(x_total)
    names(x) <- genes

    if (use_poisson_depth) {
      lam <- x * depth_scale
      x <- rpois(length(lam), lambda = lam)
      names(x) <- genes
    }

    pb_mat[, i] <- x

    if (sum(n_by_seg) > 0) P_cell[i, ] <- n_by_seg / sum(n_by_seg)
    if (sum(umi_by_seg) > 0) P_rna[i, ] <- umi_by_seg / sum(umi_by_seg)
  }

  message("Per-donor pseudo-bulk: used ", length(unique(donor_used)), " donors across ", n_mixtures, " mixtures")
  list(pb_mat = pb_mat, P_cell = P_cell, P_rna = P_rna, donor_used = donor_used)
}

# ---------
# Run pseudo-bulk test (with normalize_flag parameter)
# ---------
run_music_on_pseudobulk <- function(pb_mat, sce_obj, group.markers, clusters.type, normalize_flag = FALSE) {
  # Ensure numeric matrix genes x samples
  storage.mode(pb_mat) <- "numeric"

  # MuSiC uses common genes internally too, but we should intersect for safety
  common <- intersect(rownames(pb_mat), rownames(sce_obj))
  pb_mat2 <- pb_mat[common, , drop = FALSE]
  sce2b <- sce_obj[common, ]

  # ensure columns are character
  colData(sce2b)$segment_use <- as.character(colData(sce2b)$segment_use)
  colData(sce2b)$clusterType <- as.character(colData(sce2b)$clusterType)
  colData(sce2b)$donor_id <- as.character(colData(sce2b)$donor_id)

  # ============================================================
  # DIAGNOSTIC CHECKS: Identify potential NA sources
  # ============================================================
  message("  [DEBUG] pb_mat2 dims: ", nrow(pb_mat2), " x ", ncol(pb_mat2))
  message("  [DEBUG] pb_mat2 NAs: ", sum(is.na(pb_mat2)))
  message("  [DEBUG] pb_mat2 has zeros: ", sum(pb_mat2 == 0), " (", round(100 * sum(pb_mat2 == 0) / length(pb_mat2), 1), "%)")

  # Check if any samples have all zeros (would cause issues)
  col_sums <- colSums(pb_mat2)
  if (any(col_sums == 0)) {
    message("  [WARNING] ", sum(col_sums == 0), " pseudobulk samples have zero total counts!")
  }

  # Check SCE for issues
  sce_counts <- counts(sce2b)
  message("  [DEBUG] SCE counts NAs: ", sum(is.na(sce_counts)))

  # Check clusterType and segment_use alignment
  unique_clusterType <- unique(colData(sce2b)$clusterType)
  unique_segment <- unique(colData(sce2b)$segment_use)
  message("  [DEBUG] clusterType levels: ", paste(unique_clusterType, collapse = ", "))
  message("  [DEBUG] segment_use levels: ", paste(unique_segment, collapse = ", "))
  message("  [DEBUG] clusters.type keys: ", paste(names(clusters.type), collapse = ", "))

  # CRITICAL CHECK: clusters.type keys must match clusterType values
  missing_from_data <- setdiff(names(clusters.type), unique_clusterType)
  missing_from_spec <- setdiff(unique_clusterType, names(clusters.type))
  if (length(missing_from_data) > 0) {
    message("  [WARNING] clusters.type has keys not in data: ", paste(missing_from_data, collapse = ", "))
  }
  if (length(missing_from_spec) > 0) {
    message("  [ERROR] clusterType values missing from clusters.type: ", paste(missing_from_spec, collapse = ", "))
  }

  # Check group.markers overlap with common genes
  for (grp in names(group.markers)) {
    if (!is.null(group.markers[[grp]]) && length(group.markers[[grp]]) > 0) {
      overlap <- length(intersect(group.markers[[grp]], common))
      message("  [DEBUG] group.markers$", grp, ": ", length(group.markers[[grp]]), " genes, ", overlap, " in common")
      if (overlap == 0) {
        message("  [WARNING] group.markers$", grp, " has NO genes in common gene set!")
      }
    }
  }

  # ============================================================
  # Call MuSiC
  # ============================================================
  MuSiC::music_prop.cluster(
    bulk.mtx      = pb_mat2,
    sc.sce        = sce2b,
    group.markers = group.markers,
    groups        = "clusterType",
    clusters      = "segment_use",
    samples       = "donor_id",
    clusters.type = clusters.type,
    normalize     = normalize_flag,
    centered      = FALSE,
    verbose       = FALSE
  )
}

evaluate_pseudobulk <- function(P_true, P_hat) {
  # align columns
  common_segs <- intersect(colnames(P_true), colnames(P_hat))
  T <- P_true[, common_segs, drop = FALSE]
  H <- P_hat[, common_segs, drop = FALSE]

  # per-segment stats
  stats <- data.frame(
    segment = common_segs,
    spearman_rho = NA_real_,
    rmse = NA_real_
  )

  for (j in seq_along(common_segs)) {
    s <- common_segs[j]
    stats$spearman_rho[j] <- suppressWarnings(cor(T[, s], H[, s], method = "spearman"))
    stats$rmse[j] <- rmse(T[, s], H[, s])
  }

  list(stats = stats, T = T, H = H)
}

# ============================================================
# ACTUAL CALLS (this will run the pseudo-bulk test)
# ============================================================

message("\n================= PSEUDO-BULK TEST (EXACT) =================")

# Use improved make_pseudobulk_exact by default
EXP_USE_EXACT_PSEUDOBULK <- TRUE

if (EXP_USE_EXACT_PSEUDOBULK) {
  message("\n*** Using make_pseudobulk_exact with all fixes ***")
  message("  - P_rna_true from actual UMI contributions")
  message("  - Multinomial exact cell counts (not Dirichlet direct)")
  message("  - Independent DCT variation via two-level design")
  message("  - Multinomial downsampling for library size control")
  message("  - Per-donor sampling for MuSiC-faithful testing")

  pb <- make_pseudobulk_exact(
    sce_obj = sce2,
    segment_col = "segment_use",
    donor_col = "donor_id",
    n_mixtures = 60,
    cells_per_mixture = 2000,
    target_depth = 2000000, # Fixed depth for every mixture (2M reads)
    use_multinomial_depth = TRUE,
    depth_scale = 1.0, # irrelevant when target_depth is set
    # Independent tubule segment ranges
    dct_range = c(0.00, 0.25),
    pt_range = c(0.40, 0.75),
    tal_range = c(0.05, 0.25),
    cd_range = c(0.00, 0.15),
    other_tubule_range = c(0.00, 0.10),
    # Coarse group allocation
    tubule_total_range = c(0.50, 0.85),
    glom_total_range = c(0.05, 0.20),
    seed = 42
  )
} else if (EXP_PER_DONOR_PSEUDOBULK) {
  # Fallback: original per-donor method
  message("\n*** EXPERIMENT B: Using per-donor pseudo-bulk sampling ***")
  n_mixtures <- 60
  cells_per_mixture <- 5000
  base <- prop.table(table(colData(sce2)$segment_use))
  alpha <- as.numeric(base) * 30
  names(alpha) <- names(base)

  pb <- make_pseudobulk_per_donor(
    sce_obj = sce2,
    segment_col = "segment_use",
    donor_col = "donor_id",
    n_mixtures = n_mixtures,
    cells_per_mixture = cells_per_mixture,
    alpha = alpha,
    min_cells_per_segment = 0,
    use_poisson_depth = FALSE
  )
} else {
  # Fallback: original pooled method
  n_mixtures <- 60
  cells_per_mixture <- 5000
  base <- prop.table(table(colData(sce2)$segment_use))
  alpha <- as.numeric(base) * 30
  names(alpha) <- names(base)

  pb <- make_pseudobulk(
    sce_obj = sce2,
    n_mixtures = n_mixtures,
    cells_per_mixture = cells_per_mixture,
    alpha = alpha,
    min_cells_per_segment = 0,
    use_poisson_depth = FALSE
  )
}

# Handle different return structures
if ("P_cell_true" %in% names(pb)) {
  # New exact function returns P_cell_true and P_rna_true
  pb$P_cell <- pb$P_cell_true
  pb$P_rna <- pb$P_rna_true
}

# ============================================================
# Run MuSiC TWICE: with normalize=FALSE and normalize=TRUE
# ============================================================
message("\n================= Running MuSiC with normalize=FALSE =================")
res_F <- run_music_on_pseudobulk(pb$pb_mat, sce2, group.markers, clusters.type, normalize_flag = FALSE)

message("\n================= Running MuSiC with normalize=TRUE =================")
res_T <- run_music_on_pseudobulk(pb$pb_mat, sce2, group.markers, clusters.type, normalize_flag = TRUE)

# Evaluate both against P_rna (the correct ground truth)
eval_F <- evaluate_pseudobulk(pb$P_rna[rownames(res_F$Est.prop.weighted.cluster), ], res_F$Est.prop.weighted.cluster)
eval_T <- evaluate_pseudobulk(pb$P_rna[rownames(res_T$Est.prop.weighted.cluster), ], res_T$Est.prop.weighted.cluster)

message("\n--- normalize=FALSE: DCT recovery ---")
print(eval_F$stats[eval_F$stats$segment == "DCT", ])

message("\n--- normalize=TRUE: DCT recovery ---")
print(eval_T$stats[eval_T$stats$segment == "DCT", ])

message("\n--- Full comparison: normalize=FALSE vs normalize=TRUE ---")
compare_norm <- merge(
  eval_F$stats[, c("segment", "spearman_rho", "rmse")],
  eval_T$stats[, c("segment", "spearman_rho", "rmse")],
  by = "segment", suffixes = c("_F", "_T")
)
compare_norm$normalize_T_better <- compare_norm$spearman_rho_T > compare_norm$spearman_rho_F
print(compare_norm)

# Use the better-performing normalization for downstream
if (mean(eval_T$stats$spearman_rho, na.rm = TRUE) > mean(eval_F$stats$spearman_rho, na.rm = TRUE)) {
  message("\n*** normalize=TRUE performs better overall, using for downstream ***")
  res_pb <- res_T
  eval_rna <- eval_T
  best_normalize <- TRUE
} else {
  message("\n*** normalize=FALSE performs better overall, using for downstream ***")
  res_pb <- res_F
  eval_rna <- eval_F
  best_normalize <- FALSE
}

P_hat <- res_pb$Est.prop.weighted.cluster

# align to pseudo-bulk sample names
# MuSiC returns rows = bulk samples; our pb_mat cols are samples
# So rownames(P_hat) should match colnames(pb_mat)
stopifnot(all(rownames(P_hat) %in% colnames(pb$pb_mat)))

# ---- Evaluate against P_cell (cell fraction) ----
message("\n--- Evaluation vs P_cell (cell fraction ground truth) ---")
eval_cell <- evaluate_pseudobulk(pb$P_cell[rownames(P_hat), , drop = FALSE], P_hat)
print(eval_cell$stats[order(-eval_cell$stats$spearman_rho), ])

# Compare: if P_rna recovers well but P_cell doesn't, model is correct
message("\n--- Comparison: P_rna vs P_cell recovery ---")
comparison <- merge(
  eval_cell$stats[, c("segment", "spearman_rho")],
  eval_rna$stats[, c("segment", "spearman_rho")],
  by = "segment", suffixes = c("_cell", "_rna")
)
comparison$rna_better <- comparison$spearman_rho_rna > comparison$spearman_rho_cell
print(comparison)

if (mean(comparison$spearman_rho_rna, na.rm = TRUE) > mean(comparison$spearman_rho_cell, na.rm = TRUE) + 0.05) {
  message("\n*** MuSiC recovers P_rna better than P_cell - model is fine, evaluation was mis-specified ***")
} else {
  message("\n*** P_cell and P_rna recovery are similar ***")
}

# Save results - use P_rna as primary (correct ground truth for deconvolution)
write.csv(eval_rna$stats, file.path(outdir, "pseudobulk_recovery_stats_P_rna.csv"), row.names = FALSE)
write.csv(eval_cell$stats, file.path(outdir, "pseudobulk_recovery_stats_P_cell.csv"), row.names = FALSE)
write.csv(comparison, file.path(outdir, "pseudobulk_recovery_comparison.csv"), row.names = FALSE)
write.csv(eval_rna$T, file.path(outdir, "pseudobulk_true_P_rna.csv"), quote = FALSE)
write.csv(eval_cell$T, file.path(outdir, "pseudobulk_true_P_cell.csv"), quote = FALSE)
write.csv(eval_rna$H, file.path(outdir, "pseudobulk_est_props.csv"), quote = FALSE)
write.csv(compare_norm, file.path(outdir, "pseudobulk_normalize_comparison.csv"), row.names = FALSE)
message("Wrote pseudo-bulk outputs to: ", normalizePath(outdir))
message("Best normalize setting: ", best_normalize)

# Quick scatter plots - show P_rna recovery (the correct metric)
df_plot <- do.call(rbind, lapply(colnames(eval_rna$T), function(s) {
  data.frame(segment = s, true = eval_rna$T[, s], est = eval_rna$H[, s])
}))

p_rec <- ggplot(df_plot, aes(x = true, y = est)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~segment, scales = "free") +
  theme_bw() +
  labs(title = "Pseudo-bulk recovery: P_rna (UMI fraction) vs estimated", x = "True RNA proportion", y = "Estimated proportion")
ggsave(file.path(outdir, "pseudobulk_recovery_scatter_P_rna.png"), p_rec, width = 11, height = 7)
message("Saved pseudo-bulk recovery plot: pseudobulk_recovery_scatter_P_rna.png")

res_cluster <- MuSiC::music_prop.cluster(
  bulk.mtx      = bulk_mat2, # genes x samples
  sc.sce        = sce2,
  group.markers = group.markers,
  groups        = "clusterType",
  clusters      = "segment_use",
  samples       = "donor_id", # must exist in colData(sce2)
  clusters.type = clusters.type,
  normalize     = TRUE, # IMPORTANT (see MuSiC docs)
  centered      = FALSE,
  verbose       = TRUE
)

# Extract cell-type proportions
prop_seg <- res_cluster$Est.prop.weighted.cluster
write.csv(prop_seg, file.path(outdir, "music_cluster_proportions.csv"), quote = FALSE)
message("Wrote tree-guided cluster proportions to: ", file.path(outdir, "music_cluster_proportions.csv"))

# Compute coarse group proportions by summing segments (MuSiC doesn't return this directly)
# Safely sum segments that exist in prop_seg
safe_sum_cols <- function(mat, cols) {
  existing <- intersect(cols, colnames(mat))
  if (length(existing) == 0) {
    return(rep(0, nrow(mat)))
  }
  rowSums(mat[, existing, drop = FALSE])
}

# Build prop_group dynamically based on actual clusters.type keys
prop_group <- as.data.frame(lapply(names(clusters.type), function(grp) {
  safe_sum_cols(prop_seg, clusters.type[[grp]])
}))
colnames(prop_group) <- names(clusters.type)
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
    if (sum(ok) == 0) {
      return(character(0))
    }

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
message(
  "CPM normalized bulk matrix for validation (library sizes: ",
  round(min(lib_sizes) / 1e6, 1), "M - ", round(max(lib_sizes) / 1e6, 1), "M)"
)

# Function to compute bulk expression scores using each segment's markers (CPM-normalized)
compute_bulk_scores <- function(bulk_cpm_mat, markers_list) {
  scores <- sapply(names(markers_list), function(seg) {
    genes <- intersect(markers_list[[seg]], rownames(bulk_cpm_mat))
    if (length(genes) == 0) {
      return(rep(NA, ncol(bulk_cpm_mat)))
    }
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
    } else {
      NA
    }
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
  row.names = FALSE
)
message("\nWrote validation results to: ", file.path(outdir, "validation_segment_correlations.csv"))

# ---------------------------
# DCT-specific marker check
# ---------------------------
message("\n--- DCT Marker Validation ---")
dct_markers <- c("Slc12a3", "Pvalb", "Calb1", "Trpv5")
dct_genes <- resolve_markers_to_matrix(dct_markers, rownames(bulk_cpm))
message("DCT markers found in bulk: ", paste(dct_genes, collapse = ", "))

if (length(dct_genes) > 0) {
  dct_score <- colMeans(log1p(bulk_cpm[dct_genes, , drop = FALSE]))
  dct_cor <- cor.test(prop_seg[, "DCT"], dct_score, method = "spearman")
  message(sprintf(
    "DCT proportion vs DCT marker score: rho=%.3f, p=%.4g",
    dct_cor$estimate, dct_cor$p.value
  ))
} else {
  message("WARNING: No DCT markers found in bulk data for validation.")
}

# Summary statistics
good_segments <- validation_results$segment[!is.na(validation_results$spearman_rho) &
  validation_results$spearman_rho > 0.3]
poor_segments <- validation_results$segment[!is.na(validation_results$spearman_rho) &
  validation_results$spearman_rho < 0.1]

message("\n--- Validation Summary ---")
message(
  "Segments with good correlation (rho > 0.3): ",
  if (length(good_segments) > 0) paste(good_segments, collapse = ", ") else "None"
)
message(
  "Segments with poor correlation (rho < 0.1): ",
  if (length(poor_segments) > 0) paste(poor_segments, collapse = ", ") else "None"
)

# ---------------------------
# 5c) Coarse Group Validation
# ---------------------------
message("\n--- Coarse Group Split Sanity Check ---")

for (grp in names(clusters.type)) {
  segs <- intersect(clusters.type[[grp]], colnames(prop_seg))
  if (length(segs) == 0) next

  grp_total <- rowSums(prop_seg[, segs, drop = FALSE])
  message(sprintf("\n%s (segments: %s):", grp, paste(segs, collapse = ", ")))
  message(sprintf(
    "  Mean total: %.3f, Median: %.3f, Range: [%.4f, %.3f]",
    mean(grp_total), median(grp_total), min(grp_total), max(grp_total)
  ))

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
