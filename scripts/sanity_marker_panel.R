# =============================================================================
# Sanity Marker Panel Script
# Purpose: Validate marker expression across nephron segments to ensure correct
#          cell type labeling in the reference SCE before deconvolution.
#
# This script:
#   1) Maps marker panel genes into the gene ID space used by sce2
#   2) Computes mean logCPM and % cells expressing for each marker by segment_use
#   3) Identifies where each marker peaks to detect potential mislabeling
#   4) Optionally refines segment labels (DCT1, DCT2/CNT, CD_PC, CD_IC) using
#      module-score winner assignment for improved deconvolution
#
# Prerequisites: Run this AFTER run_deconvolution.R has created `sce2` so that
#                `sce2` is in your global environment with segment_use assigned.
# =============================================================================

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
})

stopifnot(exists("sce2"))
stopifnot("segment_use" %in% colnames(colData(sce2)))
stopifnot("counts" %in% assayNames(sce2))

# =============================================================================
# PART A: Sanity Marker Panel (means + % expressing)
# =============================================================================

# ----------------------------
# 0) Marker panel (SYMBOLS)
# ----------------------------
panel <- list(
    PT      = c("Slc34a1", "Lrp2", "Aqp1"),
    TAL     = c("Slc12a1", "Umod", "Kcnj1", "Clcnkb"),
    DCT1    = c("Slc12a3", "Pvalb"),
    DCT2CNT = c("Calb1", "Trpv5", "Slc8a1", "Klk1"),
    CD_PC   = c("Aqp2", "Avpr2", "Scnn1g"),
    CD_IC   = c("Atp6v1b1", "Slc4a1", "Foxi1")
)
all_syms <- unique(unlist(panel))

# ----------------------------
# 1) Map SYMBOL -> ENSMUSG (only if needed)
# ----------------------------
rown_is_ens <- mean(grepl("^ENSMUSG", rownames(sce2))) > 0.5

sym2ens <- function(symbols) {
    m <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = symbols,
        keytype = "SYMBOL",
        columns = c("ENSEMBL", "SYMBOL")
    )
    m <- m[!is.na(m$ENSEMBL), ]
    # keep first mapping per symbol (good enough for markers)
    m <- m[!duplicated(m$SYMBOL), ]
    setNames(m$ENSEMBL, m$SYMBOL)
}

map <- sym2ens(all_syms)

# pick marker IDs in the same space as your SCE rownames
panel_ids <- lapply(panel, function(syms) {
    if (rown_is_ens) unname(map[syms]) else syms
})
# drop NAs and keep only genes present in sce2
panel_ids <- lapply(panel_ids, function(ids) intersect(na.omit(ids), rownames(sce2)))

cat("\n=== Marker presence in reference ===\n")
for (k in names(panel_ids)) {
    cat(sprintf("%-7s: %d/%d present\n", k, length(panel_ids[[k]]), length(panel[[k]])))
    if (length(panel_ids[[k]]) == 0) cat("  !!! NONE present for this group\n")
}

# ----------------------------
# 2) Compute logCPM for scoring + mean per segment
# ----------------------------
X <- assay(sce2, "counts")
grp <- factor(colData(sce2)$segment_use)

# logCPM: X * diag(1/libsize) * 1e6, then log1p
lib <- Matrix::colSums(X)
Dinv <- Matrix::Diagonal(x = 1 / pmax(lib, 1))
logcpm <- log1p((X %*% Dinv) * 1e6)

# group indicator (cells x groups)
G <- Matrix::sparse.model.matrix(~ 0 + grp)
colnames(G) <- sub("^grp", "", colnames(G))
n_per_group <- Matrix::colSums(G)

# group mean helper: (genes x cells) %*% (cells x groups) => genes x groups
group_means <- function(M, gene_ids) {
    if (length(gene_ids) == 0) {
        return(NULL)
    }
    S <- M[gene_ids, , drop = FALSE] %*% G
    sweep(S, 2, n_per_group, "/")
}

# percent expressing helper (fraction of cells with counts>0)
group_pct_expr <- function(counts, gene_ids) {
    if (length(gene_ids) == 0) {
        return(NULL)
    }
    B <- (counts[gene_ids, , drop = FALSE] > 0)
    # convert logical sparse -> numeric sparse
    B <- Matrix::Matrix(B * 1, sparse = TRUE)
    S <- B %*% G
    sweep(S, 2, n_per_group, "/")
}

# ----------------------------
# 3) Build tables for the full panel
# ----------------------------
panel_gene_ids <- unique(unlist(panel_ids))
if (length(panel_gene_ids) == 0) stop("No marker genes found in sce2 after mapping.")

means <- group_means(logcpm, panel_gene_ids)
pct <- group_pct_expr(X, panel_gene_ids)

# display names: SYMBOL (ENSMUSG...) if Ensembl
if (rown_is_ens) {
    ens2sym <- setNames(names(map), unname(map))
    rn <- rownames(means)
    pretty <- sapply(rn, function(e) {
        s <- ens2sym[[e]]
        if (!is.null(s) && !is.na(s)) paste0(s, " (", e, ")") else e
    })
    rownames(means) <- pretty
    rownames(pct) <- pretty
}

cat("\n=== Mean logCPM by segment (marker panel) ===\n")
print(round(means, 2))

cat("\n=== % cells expressing (counts>0) by segment ===\n")
print(round(pct * 100, 1))

# ----------------------------
# 4) Quick "where does each marker peak?"
# ----------------------------
peak_seg <- apply(means, 1, function(v) colnames(means)[which.max(v)])
cat("\n=== Peak segment per marker (by mean logCPM) ===\n")
print(peak_seg)

# Helpful check: do CD principal markers peak in Other?
cdpc_rows <- if (rown_is_ens) {
    # match by symbol in display names
    grep("^Aqp2 \\(|^Avpr2 \\(|^Scnn1g \\(", rownames(means), value = TRUE)
} else {
    c("Aqp2", "Avpr2", "Scnn1g")
}

if (length(cdpc_rows) > 0) {
    cat("\n=== CD principal markers: peak segments ===\n")
    print(peak_seg[cdpc_rows])
}

# =============================================================================
# How to interpret the output (the "skeptical" read):
#
# - If Aqp2 / Avpr2 / Scnn1g peak in "Other", you likely have true collecting
#   duct principal cells hiding in "Other".
# - If Atp6v1b1 / Foxi1 / Slc4a1 peak in "Other", you likely have intercalated
#   cells hiding in "Other".
# - If your current "CD" peaks for Calb1/Trpv5/Slc8a1/Klk1, that's CNT/DCT2,
#   NOT true collecting duct.
# - If Slc12a1/Umod don't peak in your TAL_LOH, your TAL label may be off
#   (or TAL is underrepresented).
# =============================================================================


# =============================================================================
# PART B: Rename/split distal segments (practical way)
# =============================================================================
# Option 1: Simple "module-score winner" relabel inside (DCT, CD, Other)
#
# This is the fastest approach, and it's surprisingly effective.
# After running this, use `segment_refined` instead of `segment_use` in MuSiC.
# =============================================================================

# Set to TRUE to enable the relabeling (run Part A first to inspect markers)
RUN_RELABELING <- FALSE

if (RUN_RELABELING) {
    message("\n=== Running distal segment refinement ===")

    # pick the subset you want to relabel
    distal <- colData(sce2)$segment_use %in% c("DCT", "CD", "Other", "TAL_LOH")
    cells <- which(distal)

    # reuse logcpm from above (genes x cells)
    logcpm_dist <- logcpm[, cells, drop = FALSE]

    # helper: average marker expression per cell (logCPM)
    score_set <- function(ids) {
        ids <- intersect(ids, rownames(sce2))
        if (length(ids) == 0) {
            return(rep(NA_real_, ncol(logcpm_dist)))
        }
        Matrix::colMeans(logcpm_dist[ids, , drop = FALSE])
    }

    # marker IDs (from above mapping step)
    dct1_ids <- panel_ids$DCT1
    dct2_ids <- panel_ids$DCT2CNT
    cdpc_ids <- panel_ids$CD_PC
    cdic_ids <- panel_ids$CD_IC
    tal_ids <- panel_ids$TAL

    scores <- cbind(
        DCT1    = score_set(dct1_ids),
        DCT2CNT = score_set(dct2_ids),
        CD_PC   = score_set(cdpc_ids),
        CD_IC   = score_set(cdic_ids),
        TAL     = score_set(tal_ids)
    )

    # assign label = argmax score, but only if it's convincingly high
    best_lab <- colnames(scores)[max.col(scores, ties.method = "first")]
    best_val <- apply(scores, 1, max, na.rm = TRUE)

    # threshold: tune this; start conservative
    thr <- 0.3 # logCPM threshold; raise if overcalling
    new_lab <- as.character(colData(sce2)$segment_use[cells])
    new_lab[best_val >= thr] <- best_lab[best_val >= thr]

    # write back as a NEW column so you don't overwrite your old labels
    colData(sce2)$segment_refined <- as.character(colData(sce2)$segment_use)
    colData(sce2)$segment_refined[cells] <- new_lab
    colData(sce2)$segment_refined <- factor(colData(sce2)$segment_refined)

    cat("\n=== Before/after counts (distal relabeling) ===\n")
    print(table(before = colData(sce2)$segment_use[cells]))
    print(table(after = colData(sce2)$segment_refined[cells]))

    # =============================================================================
    # What this accomplishes:
    #   - If true collecting duct is hiding in "Other", it becomes CD_PC / CD_IC.
    #   - Your current "CD" that's really CNT becomes DCT2CNT.
    #   - Your "DCT" that's really DCT1 stays DCT1.
    #   - You can keep TAL separate or ignore it.
    #
    # IMPORTANT: Tell MuSiC to use `segment_refined` instead of `segment_use`
    #   Wherever your deconvolution script currently uses `segment_use`,
    #   switch the cluster column to `segment_refined`.
    # =============================================================================

    message("\nDone! Use colData(sce2)$segment_refined for MuSiC deconvolution.")
    message("Remember to update clusters.type to include: DCT1, DCT2CNT, CD_PC, CD_IC")
} else {
    message("\n=== Relabeling disabled ===")
    message("Set RUN_RELABELING <- TRUE after inspecting Part A output to enable refinement.")
}

# =============================================================================
# Summary statistics export (optional)
# =============================================================================

# Save marker expression tables for later review
outdir <- "data/processed/deconvolution/marker_sanity"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.csv(as.matrix(means), file.path(outdir, "marker_mean_logcpm_by_segment.csv"), quote = FALSE)
write.csv(as.matrix(pct * 100), file.path(outdir, "marker_pct_expressing_by_segment.csv"), quote = FALSE)
write.csv(data.frame(marker = names(peak_seg), peak_segment = peak_seg),
    file.path(outdir, "marker_peak_segments.csv"),
    row.names = FALSE
)

message("\nSaved marker sanity outputs to: ", normalizePath(outdir))
message("Files: marker_mean_logcpm_by_segment.csv, marker_pct_expressing_by_segment.csv, marker_peak_segments.csv")
