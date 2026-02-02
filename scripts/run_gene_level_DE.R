# scripts/run_gene_level_DE.R
#
# Generate gene-level differential expression tables for FLT vs GC contrasts
# within each Age×Arm stratum. These outputs can be used with
# phase5_build_silent_shifters_strict.py --gene_de
#
# Outputs: data/processed/gene_level_DE/
#   - ISS_T_YNG_FLT_vs_GC_gene_DE.tsv
#   - ISS_T_OLD_FLT_vs_GC_gene_DE.tsv
#   - LAR_YNG_FLT_vs_GC_gene_DE.tsv
#   - LAR_OLD_FLT_vs_GC_gene_DE.tsv

suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(data.table)
})

# Check for ashr (optional, for LFC shrinkage)
HAS_ASHR <- requireNamespace("ashr", quietly = TRUE)
if (HAS_ASHR) {
    cat("ashr package found - will use LFC shrinkage\n")
} else {
    cat("ashr package not found - skipping LFC shrinkage\n")
}

numify <- function(x) as.numeric(gsub("[^0-9eE\\.+-]", "", as.character(x)))
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

# 0) Load input data
cat("Loading bulk counts and metadata...\n")

bulk_dt <- fread("data/processed/aligned_outputs/rsem_rRNArm_raw_counts.csv")
gene_col <- names(bulk_dt)[1]
bulk_genes <- bulk_dt[[gene_col]]
bulk_counts <- as.matrix(bulk_dt[, -1, with = FALSE])
rownames(bulk_counts) <- bulk_genes
storage.mode(bulk_counts) <- "numeric"

# Strip Ensembl versions and collapse duplicates
rownames(bulk_counts) <- strip_ensembl_version(rownames(bulk_counts))
if (anyDuplicated(rownames(bulk_counts)) > 0) {
    cat("Collapsing duplicate gene IDs by sum...\n")
    bulk_counts <- rowsum(bulk_counts, group = rownames(bulk_counts), reorder = FALSE)
}

cat("Bulk counts:", nrow(bulk_counts), "genes x", ncol(bulk_counts), "samples\n")

# Load metadata
meta <- as.data.frame(fread("data/processed/aligned_outputs/metadata_aligned.tsv", header = TRUE))

# Ensure output directory exists
out_dir <- "data/processed/gene_level_DE"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1) Align meta to count columns
clean_id <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    x <- gsub("\\.(bam|sam|fastq|fq|gz|txt|csv)$", "", x, ignore.case = TRUE)
    x <- basename(x)
    x <- gsub("[\\./\\-]", "_", x)
    x <- gsub("_+", "_", x)
    x
}

count_names_raw <- colnames(bulk_counts)
meta_names_raw <- meta[["Sample Name"]]
count_names <- clean_id(count_names_raw)
meta_names <- clean_id(meta_names_raw)

# Guardrail: check for duplicates after cleaning
if (anyDuplicated(count_names)) {
    stop("Duplicate cleaned count column names; clean_id too aggressive.")
}
if (anyDuplicated(meta_names)) {
    stop("Duplicate cleaned meta sample names; clean_id too aggressive.")
}

meta_map <- setNames(seq_along(meta_names), meta_names)
idx <- meta_map[count_names]

if (any(is.na(idx))) {
    stop("Failed to align metadata to counts")
}

meta2 <- meta[idx, , drop = FALSE]
meta2[["Sample Name (raw_counts_colname)"]] <- count_names_raw
rownames(meta2) <- count_names_raw

stopifnot(nrow(meta2) == ncol(bulk_counts))
cat("Aligned", nrow(meta2), "samples\n")

# 2) Normalize factor labels (FIXED: don't collapse groups!)
normalize_age <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x <- gsub("^YOUNG$|^YNG$", "YNG", x)
    x <- gsub("^OLD$", "OLD", x)
    factor(x, levels = c("YNG", "OLD"))
}

normalize_arm <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x <- gsub("^ISS$|^ISST$|^ISS_T$|^ISS.T$", "ISS-T", x)
    x <- gsub("^LAR_T$|^LAR-T$|^LAR.T$|^LAR$", "LAR", x)
    factor(x, levels = c("ISS-T", "LAR"))
}

normalize_env <- function(x) {
    x <- toupper(trimws(as.character(x)))

    # Standardize labels WITHOUT collapsing distinct groups
    x <- gsub("^GROUND CONTROL$|^GROUND_CONTROL$", "GC", x)
    x <- gsub("^HABITAT GROUND CONTROL$|^HABITAT_GROUND_CONTROL$|^HGC/GC$", "HGC", x)
    x <- gsub("^VIVARIUM$", "VIV", x)
    x <- gsub("^VGC$|^VIV/VGC$", "VIV", x) # VGC -> VIV if those exist
    x <- gsub("^FLIGHT$|^ISS$", "FLT", x)

    # Keep FLT, GC, HGC, VIV as distinct levels
    factor(x, levels = c("GC", "HGC", "VIV", "FLT"))
}

meta2$Age <- normalize_age(meta2$Age)
meta2$Arm <- normalize_arm(meta2$Arm)
meta2$EnvGroup <- normalize_env(meta2$EnvGroup)

cat("\nFactor levels after normalization:\n")
cat("  Age:", paste(levels(meta2$Age), collapse = ", "), "\n")
cat("  Arm:", paste(levels(meta2$Arm), collapse = ", "), "\n")
cat("  EnvGroup:", paste(levels(meta2$EnvGroup), collapse = ", "), "\n")

# Check for any NAs introduced by factor levels
if (any(is.na(meta2$Age))) {
    cat("WARNING: Some Age values didn't match expected levels\n")
    print(table(meta$Age, useNA = "ifany"))
}
if (any(is.na(meta2$Arm))) {
    cat("WARNING: Some Arm values didn't match expected levels\n")
    print(table(meta$Arm, useNA = "ifany"))
}
if (any(is.na(meta2$EnvGroup))) {
    cat("WARNING: Some EnvGroup values didn't match expected levels:\n")
    print(table(meta$EnvGroup, useNA = "ifany"))
}

# 3) Create DESeqDataSet with cell-means model
cat("\nCreating DESeqDataSet...\n")

# Create cell factor (Age:Arm:EnvGroup)
meta2$cell <- factor(paste(meta2$Age, meta2$Arm, meta2$EnvGroup, sep = "_"))
cat("Cell groups:", paste(levels(meta2$cell), collapse = ", "), "\n")

# Count samples per cell (IMPORTANT: verify FLT vs GC are distinct!)
cell_counts <- table(meta2$cell)
cat("\nSamples per cell:\n")
print(cell_counts)

# Check that we have the expected cells for contrasts
expected_cells <- c(
    "YNG_ISS-T_FLT", "YNG_ISS-T_GC",
    "OLD_ISS-T_FLT", "OLD_ISS-T_GC",
    "YNG_LAR_FLT", "YNG_LAR_GC",
    "OLD_LAR_FLT", "OLD_LAR_GC"
)
missing_cells <- setdiff(expected_cells, names(cell_counts))
if (length(missing_cells) > 0) {
    cat("\nWARNING: Missing expected cells:", paste(missing_cells, collapse = ", "), "\n")
}

dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(bulk_counts)),
    colData = meta2,
    design = ~cell
)

# CPM filter
cpm_mat <- edgeR::cpm(counts(dds))
keep <- rowMeans(cpm_mat >= 1) >= 0.20
cat("\nKeeping", sum(keep), "of", length(keep), "genes after CPM filter\n")
dds <- dds[keep, ]

# Run DESeq with outlier replacement disabled (important for small n)
cat("Running DESeq2 (with minReplicatesForReplace=Inf for small n)...\n")
dds <- DESeq(dds, minReplicatesForReplace = Inf)

# 4) Extract contrasts for each Age×Arm: FLT vs GC
cat("\nExtracting FLT vs GC contrasts...\n")

contrasts <- list(
    list(name = "ISS_T_YNG_FLT_vs_GC", numerator = "YNG_ISS-T_FLT", denominator = "YNG_ISS-T_GC"),
    list(name = "ISS_T_OLD_FLT_vs_GC", numerator = "OLD_ISS-T_FLT", denominator = "OLD_ISS-T_GC"),
    list(name = "LAR_YNG_FLT_vs_GC", numerator = "YNG_LAR_FLT", denominator = "YNG_LAR_GC"),
    list(name = "LAR_OLD_FLT_vs_GC", numerator = "OLD_LAR_FLT", denominator = "OLD_LAR_GC")
)

available_cells <- levels(meta2$cell)
cat("Available cell levels:", paste(available_cells, collapse = ", "), "\n")

for (ctr in contrasts) {
    cat("\n--- Contrast:", ctr$name, "---\n")
    cat("  Numerator:", ctr$numerator, "\n")
    cat("  Denominator:", ctr$denominator, "\n")

    # Check if both levels exist
    if (!(ctr$numerator %in% available_cells)) {
        cat("  [SKIP] Numerator cell not found in data\n")
        next
    }
    if (!(ctr$denominator %in% available_cells)) {
        cat("  [SKIP] Denominator cell not found in data\n")
        next
    }

    # Print sample counts for this contrast
    n_num <- sum(meta2$cell == ctr$numerator)
    n_den <- sum(meta2$cell == ctr$denominator)
    cat("  Sample counts: FLT =", n_num, ", GC =", n_den, "\n")

    # Extract results
    res <- results(dds,
        contrast = c("cell", ctr$numerator, ctr$denominator),
        alpha = 0.05
    )

    # Apply LFC shrinkage if ashr is available (conservative log2FC for silent shifters)
    if (HAS_ASHR) {
        cat("  Applying LFC shrinkage (ashr)...\n")
        res <- lfcShrink(dds,
            contrast = c("cell", ctr$numerator, ctr$denominator),
            res = res,
            type = "ashr"
        )
    }

    # Convert to data.frame
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    # Rename columns for compatibility with silent shifters script
    res_df <- res_df[, c("gene", "log2FoldChange", "padj")]
    colnames(res_df) <- c("gene", "log2FC", "FDR")

    # Remove NA FDR (genes with too low counts)
    res_df <- res_df[!is.na(res_df$FDR), ]

    # Sort by FDR
    res_df <- res_df[order(res_df$FDR), ]

    # Write output
    out_file <- file.path(out_dir, paste0(ctr$name, "_gene_DE.tsv"))
    write.table(res_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

    cat("  Wrote:", out_file, "\n")
    cat("  Genes:", nrow(res_df), "\n")
    cat("  Significant (FDR < 0.05):", sum(res_df$FDR < 0.05), "\n")
    cat("  |log2FC| > 1:", sum(abs(res_df$log2FC) > 1), "\n")
    cat(
        "  Silent shifter candidates (|log2FC| < 0.3 & FDR > 0.2):",
        sum(abs(res_df$log2FC) < 0.3 & res_df$FDR > 0.2), "\n"
    )
}

cat("\n=== Gene-Level DE Complete ===\n")
cat("Outputs saved to:", out_dir, "\n")
