# Quick script to check DCT marker overlap with TAL/CD
# Run this in R after loading your deconvolution environment

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(SingleCellExperiment)
})

# Load the reference data
raw_dir <- "data/processed/deconvolution/raw_ref"
mtx_path <- file.path(raw_dir, "matrix.mtx")
obs_path <- file.path(raw_dir, "barcodes.tsv")
var_path <- file.path(raw_dir, "features.tsv")

raw_mat <- readMM(mtx_path)
obs <- data.frame(fread(obs_path, header = TRUE), row.names = 1, check.names = FALSE)
var <- data.frame(fread(var_path, header = TRUE), row.names = 1, check.names = FALSE)

if (ncol(raw_mat) == nrow(var)) raw_mat <- t(raw_mat)
colnames(raw_mat) <- rownames(obs)
rownames(raw_mat) <- rownames(var)

sce <- SingleCellExperiment(assays = list(counts = raw_mat), colData = obs)

# Segment mapping function (simplified for tubules)
segment_from_ct <- function(x) {
    x2 <- tolower(trimws(x))
    if (grepl("distal.*convoluted|distal.*tubule|\\bdct\\b", x2)) {
        return("DCT")
    }
    if (grepl("thick.*ascending|loop\\s+of\\s+henle|henle|\\btal\\b", x2)) {
        return("TAL_LOH")
    }
    if (grepl("collecting\\s+duct|principal\\s+cell|intercalated", x2)) {
        return("CD")
    }
    if (grepl("proximal.*tubule|proximal\\s+tube|\\bpt\\b", x2)) {
        return("PT")
    }
    "Other"
}

ct <- as.character(colData(sce)[["free_annotation"]])
ct[is.na(ct) | ct == "" | tolower(ct) %in% c("nan", "na", "unknown")] <- NA
keep <- !is.na(ct)
sce <- sce[, keep]
seg <- vapply(ct[keep], segment_from_ct, character(1))
colData(sce)[["segment_use"]] <- factor(seg)

# Get top markers per segment
get_top_markers <- function(sce_obj, seg_name, top_n = 30, min_lfc = 0.5) {
    X <- assay(sce_obj, "counts")
    cl <- as.character(colData(sce_obj)[["segment_use"]])
    segs <- unique(cl[cl != "Other"])

    means <- sapply(segs, function(s) Matrix::rowMeans(X[, cl == s, drop = FALSE]))

    mu_s <- means[, seg_name]
    mu_other <- rowMeans(means[, setdiff(segs, seg_name), drop = FALSE])
    lfc <- log2((mu_s + 1) / (mu_other + 1))

    ok <- is.finite(lfc) & (mu_s > 0.05) & (lfc > min_lfc)
    names(sort(lfc[ok], decreasing = TRUE))[seq_len(min(top_n, sum(ok)))]
}

cat("\n=== Getting Top 30 Markers per Tubular Segment ===\n")
dct_markers <- get_top_markers(sce, "DCT", top_n = 30)
tal_markers <- get_top_markers(sce, "TAL_LOH", top_n = 30)
cd_markers <- get_top_markers(sce, "CD", top_n = 30)
pt_markers <- get_top_markers(sce, "PT", top_n = 30)

cat("\n--- DCT Top 30 Markers ---\n")
print(dct_markers)

cat("\n--- TAL_LOH Top 30 Markers ---\n")
print(tal_markers)

cat("\n--- CD Top 30 Markers ---\n")
print(cd_markers)

cat("\n=== OVERLAP ANALYSIS ===\n")
cat("DCT ∩ TAL_LOH:", length(intersect(dct_markers, tal_markers)), "genes\n")
cat("  Genes:", paste(intersect(dct_markers, tal_markers), collapse = ", "), "\n")

cat("\nDCT ∩ CD:", length(intersect(dct_markers, cd_markers)), "genes\n")
cat("  Genes:", paste(intersect(dct_markers, cd_markers), collapse = ", "), "\n")

cat("\nTAL_LOH ∩ CD:", length(intersect(tal_markers, cd_markers)), "genes\n")
cat("  Genes:", paste(intersect(tal_markers, cd_markers), collapse = ", "), "\n")

cat("\n=== DCT MARKER SPECIFICITY ===\n")
dct_unique <- setdiff(dct_markers, union(tal_markers, cd_markers))
cat("DCT-unique markers (not in TAL or CD top-30):", length(dct_unique), "\n")
print(dct_unique)

# Also check known DCT markers
known_dct <- c("Slc12a3", "Pvalb", "Calb1", "Trpv5")
cat("\n=== Known DCT Markers Status ===\n")
for (g in known_dct) {
    in_dct <- g %in% dct_markers
    in_tal <- g %in% tal_markers
    in_cd <- g %in% cd_markers
    cat(sprintf("%s: in DCT top-30=%s, in TAL=%s, in CD=%s\n", g, in_dct, in_tal, in_cd))
}

# Expression heatmap data (for visualization)
cat("\n=== Mean Expression of DCT Markers Across Segments ===\n")
X <- assay(sce, "counts")
cl <- as.character(colData(sce)[["segment_use"]])
segs <- c("DCT", "TAL_LOH", "CD", "PT")
means <- sapply(segs, function(s) Matrix::rowMeans(X[, cl == s, drop = FALSE]))

# Show mean expression of DCT markers across segments
dct_expr <- means[dct_markers, ]
print(round(dct_expr, 2))

cat("\n=== Fold-change of DCT markers (DCT vs each other segment) ===\n")
for (other in c("TAL_LOH", "CD", "PT")) {
    fc <- (dct_expr[, "DCT"] + 1) / (dct_expr[, other] + 1)
    cat(sprintf("\nDCT vs %s (median FC = %.2fx):\n", other, median(fc)))
    print(round(sort(fc, decreasing = TRUE), 2))
}
