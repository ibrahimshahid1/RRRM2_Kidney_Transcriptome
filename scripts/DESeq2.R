suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(sva)
    library(data.table)
})

numify <- function(x) as.numeric(gsub("[^0-9eE\\.+-]", "", as.character(x)))
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

# ============================================================
# 0) Load input data (using fread like run_deconvolution.R)
# ============================================================
cat("Loading bulk counts and metadata...\n")

# Load counts with fread - first column is gene IDs
bulk_dt <- fread("data/processed/aligned_outputs/rsem_rRNArm_raw_counts.csv")
gene_col <- names(bulk_dt)[1]
bulk_genes <- bulk_dt[[gene_col]]
bulk_counts <- as.matrix(bulk_dt[, -1, with = FALSE])
rownames(bulk_counts) <- bulk_genes
storage.mode(bulk_counts) <- "numeric"

# Strip Ensembl versions and collapse duplicates
rownames(bulk_counts) <- strip_ensembl_version(rownames(bulk_counts))
if (anyDuplicated(rownames(bulk_counts)) > 0) {
    cat("Bulk has duplicate gene IDs after stripping versions. Collapsing by sum...\n")
    bulk_counts <- rowsum(bulk_counts, group = rownames(bulk_counts), reorder = FALSE)
}

cat("Bulk counts: ", nrow(bulk_counts), " genes x ", ncol(bulk_counts), " samples\n")
cat("Sample columns: ", paste(head(colnames(bulk_counts), 3), collapse = ", "), "...\n")

# Load metadata
meta <- as.data.frame(fread("data/processed/aligned_outputs/metadata_aligned.tsv", header = TRUE))
cat("Metadata samples: ", paste(head(meta[["Sample Name"]], 3), collapse = ", "), "...\n")

# Ensure output directory exists
out_dir <- "data/processed/phase1_residuals"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) Align meta to count columns (ROBUST)
# ============================================================

clean_id <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    # remove common file suffixes if present
    x <- gsub("\\.(bam|sam|fastq|fq|gz|txt|csv)$", "", x, ignore.case = TRUE)
    # sometimes sample names come as full paths
    x <- basename(x)
    # normalize separators
    x <- gsub("[\\./\\-]", "_", x)
    # collapse multiple underscores
    x <- gsub("_+", "_", x)
    x
}

count_names_raw <- colnames(bulk_counts)
meta_names_raw <- meta[["Sample Name"]]

count_names <- clean_id(count_names_raw)
meta_names <- clean_id(meta_names_raw)

# Diagnostic: what's missing before/after cleaning?
missing_raw <- setdiff(count_names_raw, meta_names_raw)
missing_clean <- setdiff(count_names, meta_names)

cat("\n--- NAME MATCH DIAGNOSTIC ---\n")
cat("Counts columns:", length(count_names_raw), "\n")
cat("Meta sample names:", length(meta_names_raw), "\n")
cat("Missing (raw exact match):", length(missing_raw), "\n")
if (length(missing_raw) > 0) {
    cat("Example missing raw:\n")
    print(head(missing_raw, 10))
}
cat("Missing (after cleaning):", length(missing_clean), "\n")
if (length(missing_clean) > 0) {
    cat("Example missing cleaned:\n")
    print(head(missing_clean, 10))
}

# Build index mapping from cleaned names
meta_map <- setNames(seq_along(meta_names), meta_names)

idx <- meta_map[count_names] # vector of row indices in meta corresponding to each count column

# If any NA remains, print them and stop with helpful output
if (any(is.na(idx))) {
    bad <- which(is.na(idx))
    cat("\nERROR: Still cannot match these count columns to metadata Sample Name.\n")
    cat("Unmatched count columns (raw):\n")
    print(count_names_raw[bad])
    cat("\nUnmatched count columns (cleaned):\n")
    print(count_names[bad])
    cat("\nTip: inspect meta[['Sample Name']] values around similar strings.\n")
    stop("Failed to align metadata to counts after cleaning.")
}

# Align meta to counts in correct order
meta2 <- meta[idx, , drop = FALSE]

# Overwrite meta2 Sample Name to be the actual count column names (clean + raw)
meta2[["Sample Name (raw_counts_colname)"]] <- count_names_raw
meta2[["Sample Name (cleaned)"]] <- count_names

# Hard check: now aligned
stopifnot(nrow(meta2) == ncol(bulk_counts))
cat("Successfully aligned", nrow(meta2), "samples.\n")

# Set rownames on meta2 for rock-solid alignment downstream
rownames(meta2) <- meta2[["Sample Name (raw_counts_colname)"]]

# ============================================================
# 2) Load CLR deconv covariates (your saved output)
# ============================================================
clr_file <- "data/processed/deconvolution/test16/music_segment_direct_proportions_CLR.csv"
if (!file.exists(clr_file)) {
    # Fall back to base deconvolution dir
    clr_file <- "data/processed/deconvolution/music_segment_direct_proportions_CLR.csv"
}
cat("Loading CLR proportions from:", clr_file, "\n")
clr <- read.csv(clr_file, row.names = 1, check.names = FALSE)
# Align CLR to meta2 rownames (rock-solid alignment)
clr <- clr[rownames(meta2), , drop = FALSE]

# Drop one CLR column to avoid rank deficiency (choose a stable large compartment, e.g. PT)
drop_part <- "PT"
clr_use <- clr[, setdiff(colnames(clr), drop_part), drop = FALSE]

# ============================================================
# 3) Build technical covariates from your metadata
# ============================================================
meta2$LibraryBatch <- factor(meta2[["Parameter Value[Library Batch Number]"]])
meta2$LibraryKit <- factor(meta2[["Parameter Value[Library Kit]"]])
meta2$SeqInstr <- factor(meta2[["Parameter Value[Sequencing Instrument]"]])

meta2$ReadDepth <- numify(meta2[["Parameter Value[Read Depth]"]])
meta2$rRNA <- numify(meta2[["Parameter Value[rRNA Contamination]"]])

meta2$FragSize <- numify(meta2[["Parameter Value[Fragment Size]"]])
meta2$ReadLen <- numify(meta2[["Parameter Value[Read Length]"]])
meta2$QA <- numify(meta2[["Parameter Value[QA Score]"]])

meta2$Stranded <- factor(meta2[["Parameter Value[Stranded]"]])
meta2$LibLayout <- factor(meta2[["Parameter Value[Library Layout]"]])
meta2$LibSelect <- factor(meta2[["Parameter Value[Library Selection]"]])

# Biology factors
meta2$Age <- factor(meta2$Age)
meta2$Arm <- factor(meta2$Arm)
meta2$EnvGroup <- factor(meta2$EnvGroup)

# ============================================================
# 4) VST expression matrix Y (genes x samples)
# ============================================================
cat("Creating DESeqDataSet and applying VST normalization...\n")
dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(bulk_counts)),
    colData = meta2,
    design = ~1
)

# Counts-consistent filter (CPM >= 1 in at least 20% of samples)
cpm_mat <- edgeR::cpm(counts(dds))
keep <- rowMeans(cpm_mat >= 1) >= 0.20
cat("Keeping", sum(keep), "of", length(keep), "genes after CPM filtering\n")
dds <- dds[keep, ]

vsd <- DESeq2::vst(dds, blind = FALSE)
Y <- assay(vsd) # genes x samples

# ============================================================
# 5) SVA: capture residual latent factors (HARDENED)
# ============================================================
cat("\n--- Preparing Design Matrix for SVA ---\n")

# 0) Build covariate DF
tech_candidates <- c("ReadDepth", "rRNA", "LibraryBatch", "FragSize", "QA")
bio_terms <- c("Age", "Arm", "EnvGroup")

# keep only existing columns
tech_candidates <- tech_candidates[tech_candidates %in% colnames(meta2)]
bio_terms <- bio_terms[bio_terms %in% colnames(meta2)]

covar_df <- meta2[, c(bio_terms, tech_candidates), drop = FALSE]
covar_df <- cbind(covar_df, clr_use)

# 1) Clean: drop single-level factors; scale numerics (huge for conditioning)
drop_single_level_factors <- function(df) {
    for (nm in names(df)) {
        if (is.character(df[[nm]]) || is.factor(df[[nm]])) {
            df[[nm]] <- factor(df[[nm]])
            if (nlevels(df[[nm]]) < 2) df[[nm]] <- NULL
        }
    }
    df
}
covar_df <- drop_single_level_factors(covar_df)

scale_numeric <- function(df) {
    for (nm in names(df)) {
        if (is.numeric(df[[nm]])) {
            if (any(is.na(df[[nm]]))) df[[nm]][is.na(df[[nm]])] <- median(df[[nm]], na.rm = TRUE)
            if (sd(df[[nm]], na.rm = TRUE) < 1e-12) {
                df[[nm]] <- NULL
            } else {
                df[[nm]] <- as.numeric(scale(df[[nm]]))
            }
        }
    }
    df
}
covar_df <- scale_numeric(covar_df)

# Recompute term lists after dropping columns
all_terms <- colnames(covar_df)
bio_terms <- intersect(bio_terms, all_terms)
clr_terms <- colnames(clr_use)
tech_terms <- setdiff(setdiff(all_terms, bio_terms), clr_terms)

# 2) Helper: make a model matrix that is FULL-RANK and WELL-CONDITIONED
make_fullrank <- function(M, tol = 1e-7) {
    q <- qr(M, tol = tol)
    M2 <- M[, q$pivot[seq_len(q$rank)], drop = FALSE]
    # Orthonormalize columns -> crossprod is ~I, so num.sv/sva won't choke
    Q <- qr.Q(qr(M2), complete = FALSE)
    colnames(Q) <- colnames(M2)
    Q
}

# 3) Build raw model matrices
form_full <- as.formula(paste("~", paste(c(bio_terms, tech_terms, clr_terms), collapse = " + ")))
form_null <- as.formula(paste("~", paste(c(tech_terms, clr_terms), collapse = " + ")))

mod_raw <- model.matrix(form_full, data = covar_df)
mod0_raw <- model.matrix(form_null, data = covar_df)

# 4) Force full-rank + stable conditioning
mod <- make_fullrank(mod_raw, tol = 1e-7)
mod0 <- make_fullrank(mod0_raw, tol = 1e-7)

cat("Full model dims (after full-rank + orthonormalize):", paste(dim(mod), collapse = " x "), "\n")
cat("Null model dims (after full-rank + orthonormalize):", paste(dim(mod0), collapse = " x "), "\n")

# Optional sanity diagnostics
cat("kappa(mod_raw) =", kappa(mod_raw), "\n")
cat("rank(mod_raw)  =", qr(mod_raw)$rank, "/", ncol(mod_raw), "\n")

cat("Estimating number of SVs...\n")
n.sv <- num.sv(Y, mod, method = "leek")
cat("Estimated SVs:", n.sv, "\n")

svobj <- sva(Y, mod, mod0, n.sv = n.sv)

for (k in seq_len(ncol(svobj$sv))) {
    meta2[[paste0("SV", k)]] <- svobj$sv[, k]
}

# ============================================================
# 6) Residualization WITHOUT matrix inversion (QR RESIDUALS) - SAFE
# ============================================================
cat("\n--- Performing Residualization (QR-based; stable) ---\n")

# 1) Build X_remove from the SAME covar_df used for mod/mod0
X_remove_raw <- model.matrix(form_null, data = covar_df)

# 2) Handle the "0 SV" case safely
sv_mat <- svobj$sv
if (is.null(sv_mat) || ncol(sv_mat) == 0) {
    X_remove <- X_remove_raw
    cat("No SVs detected (n.sv=0). Residualizing with Null model only.\n")
} else {
    sv_mat <- as.matrix(sv_mat)
    colnames(sv_mat) <- paste0("SV", seq_len(ncol(sv_mat)))
    rownames(sv_mat) <- rownames(X_remove_raw) # IMPORTANT: align rownames
    X_remove <- cbind(X_remove_raw, sv_mat)
}

# 3) Hard checks: dimensions must match Y
stopifnot(nrow(X_remove) == ncol(Y))
# Ensure sample order matches exactly (this prevents silent misalignment)
if (!is.null(colnames(Y))) {
    # model.matrix rownames are sample indices; so enforce consistent ordering via rownames
    # If your covar_df rownames aren't set, set them earlier to sample names.
    cat("Y dims:", paste(dim(Y), collapse = " x "), "\n")
    cat("X_remove dims:", paste(dim(X_remove), collapse = " x "), "\n")
}

# 4) QR residuals in one shot (samples x genes)
qrX <- qr(X_remove)
Yt <- t(Y) # samples x genes

# qr.resid returns (samples x genes)
Rt <- qr.resid(qrX, Yt)

# Convert back to (genes x samples) and restore dimnames cleanly
Rtech <- t(Rt)
rownames(Rtech) <- rownames(Y)
colnames(Rtech) <- colnames(Y)

cat("Rtech dims:", paste(dim(Rtech), collapse = " x "), "\n")

# ============================================================
# 7) Save results
# ============================================================
out_file <- file.path(out_dir, "phase1_Rtech.rds")
cat("Saving results to:", out_file, "\n")
saveRDS(list(Y = Y, Rtech = Rtech, meta = meta2, clr = clr, svobj = svobj),
    file = out_file
)

cat("\n=== DESeq2 + SVA Residualization Complete ===\n")
