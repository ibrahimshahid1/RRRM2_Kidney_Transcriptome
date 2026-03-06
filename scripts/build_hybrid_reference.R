#!/usr/bin/env Rscript
# =============================================================================
# build_hybrid_reference.R
#
# Build a hybrid single-cell reference atlas for MuSiC deconvolution by
# combining:
#   - TMS (Tabula Muris Senis): background kidney cell types (PT, Podocyte,
#     Endothelial, Mesangial, Fibroblast, Immune, CD)
#   - Chen et al. 2021 (GSE150338): high-resolution distal nephron types
#     (DCT1, DCT2, CNT, CTAL)
#
# Cluster-to-celltype mapping for Chen (from marker expression analysis):
#   DCT1  = clusters 0, 1, 2, 8, 12  (Slc12a3+, Pvalb high)
#   DCT2  = cluster 3                 (Slc12a3+, Pvalb low, Calb1+)
#   CNT   = cluster 6                 (Slc12a3+, Calb1++, Scnn1g+, Pvalb-)
#   CTAL  = clusters 4, 5, 13         (Slc12a1+, Umod+)
#   DROP  = clusters 7, 9, 10, 11     (contaminants: stromal, CD-PC, PT, IC)
#
# Output: data/processed/deconvolution/hybrid_ref/
#   matrix.mtx, barcodes.tsv, features.tsv
# =============================================================================

suppressPackageStartupMessages({
    library(Matrix)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
})

REPO <- normalizePath(".")
TMS_DIR <- file.path(REPO, "data/processed/deconvolution/raw_ref")
CHEN_RDATA <- file.path(
    REPO, "data/external/single_cell_atlases/chen_atlas",
    "GSE150338_Seurat_IntegrateData.RData"
)
OUT_DIR <- file.path(REPO, "data/processed/deconvolution/hybrid_ref")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Chen atlas ──────────────────────────────────────────────────────
message(" Loading Chen atlas…")
load(CHEN_RDATA)
chen_obj <- merge.data.integrated
chen_counts <- chen_obj@assays$RNA@counts # dgCMatrix: 17655 genes × 9099 cells
chen_meta <- chen_obj@meta.data
chen_clusters <- chen_meta$seurat_clusters

message("  Chen: ", ncol(chen_counts), " cells × ", nrow(chen_counts), " genes")

# Assign cell types based on clusters
chen_cluster_map <- c(
    "0"  = "DCT1",
    "1"  = "DCT1",
    "2"  = "DCT1",
    "3"  = "DCT2",
    "4"  = "TAL_LOH",
    "5"  = "TAL_LOH",
    "6"  = "CNT",
    "7"  = NA, # stromal/unresolved – drop
    "8"  = "DCT1",
    "9"  = NA, # CD-PC contaminant – drop
    "10" = NA, # PT contaminant – drop
    "11" = NA, # IC contaminant – drop
    "12" = "DCT1",
    "13" = "TAL_LOH"
)

chen_meta$cell_type <- chen_cluster_map[as.character(chen_clusters)]

# Drop contaminant clusters
keep <- !is.na(chen_meta$cell_type)
chen_counts <- chen_counts[, keep]
chen_meta <- chen_meta[keep, ]
message("  After filtering: ", ncol(chen_counts), " cells retained")
message("  Cell types: ")
print(table(chen_meta$cell_type))

# Assign donor_id from Batch_10X (rep1/rep2/rep3 → Chen_rep1/Chen_rep2/Chen_rep3)
chen_meta$donor_id <- paste0("Chen_", chen_meta$Batch_10X)
message("  Donors: ", paste(unique(chen_meta$donor_id), collapse = ", "))

# ── 2. Convert Chen gene symbols to Ensembl IDs ─────────────────────────────
message("\n▶ Converting gene symbols → Ensembl IDs…")
chen_symbols <- rownames(chen_counts)

ens_ids <- mapIds(org.Mm.eg.db,
    keys = chen_symbols,
    keytype = "SYMBOL",
    column = "ENSEMBL",
    multiVals = "first"
)

mapped <- !is.na(ens_ids)
message(
    "  Mapped: ", sum(mapped), " / ", length(chen_symbols),
    " symbols  (", round(100 * mean(mapped), 1), "%)"
)

chen_counts <- chen_counts[mapped, ]
chen_ensembl <- ens_ids[mapped]
rownames(chen_counts) <- chen_ensembl

# Handle duplicate Ensembl IDs (multiple symbols mapping to same Ensembl)
dup_ens <- duplicated(rownames(chen_counts))
if (any(dup_ens)) {
    message("  Removing ", sum(dup_ens), " duplicate Ensembl IDs")
    chen_counts <- chen_counts[!dup_ens, ]
}
message("  Chen final: ", nrow(chen_counts), " genes × ", ncol(chen_counts), " cells")

# ── 3. Load TMS reference ───────────────────────────────────────────────────
message("\n▶ Loading TMS reference…")
tms_barcodes <- read.delim(file.path(TMS_DIR, "barcodes.tsv"),
    stringsAsFactors = FALSE
)
tms_features <- read.delim(file.path(TMS_DIR, "features.tsv"),
    stringsAsFactors = FALSE, row.names = 1
)
tms_counts <- readMM(file.path(TMS_DIR, "matrix.mtx"))
tms_counts <- as(tms_counts, "CsparseMatrix")

# TMS MTX is stored as cells × genes; transpose to genes × cells
if (nrow(tms_counts) == nrow(tms_barcodes) &&
    ncol(tms_counts) == nrow(tms_features)) {
    message("  (transposing: cells×genes → genes×cells)")
    tms_counts <- t(tms_counts)
}

# Set dimension names
rownames(tms_counts) <- rownames(tms_features)
colnames(tms_counts) <- tms_barcodes$index

message("  TMS: ", ncol(tms_counts), " cells × ", nrow(tms_counts), " genes")

# ── 4. Remove DCT and TAL_LOH from TMS ──────────────────────────────────────
# We need to identify which TMS cells are DCT or TAL_LOH.
# Use the same segment_from_ct mapping logic from deconvolution.R
segment_from_ct <- function(ct) {
    ct <- tolower(ct)
    seg <- rep(NA_character_, length(ct))
    seg[grepl("proximal", ct)] <- "PT"
    seg[grepl("podocyte", ct)] <- "Podocyte"
    seg[grepl("endotheli", ct)] <- "Endothelial"
    seg[grepl("mesangi", ct)] <- "Mesangial"
    seg[grepl("fibroblast", ct)] <- "Fibroblast"
    seg[grepl("loop of henle|ascending limb|thick ascending", ct)] <- "TAL_LOH"
    seg[grepl("distal convoluted", ct)] <- "DCT"
    seg[grepl("collecting duct|principal cell|intercalat", ct)] <- "CD"
    seg[grepl("macrophage|monocyte|leukocyte|lymphocyte|t cell|b cell|nk cell|dendritic|immune", ct)] <- "Immune"
    seg
}

tms_segments <- segment_from_ct(tms_barcodes$free_annotation)
tms_barcodes$segment <- tms_segments

message("  TMS segment distribution before filter:")
print(table(tms_segments, useNA = "ifany"))

# Keep everything EXCEPT DCT and TAL_LOH from TMS
tms_keep <- !(tms_segments %in% c("DCT", "TAL_LOH")) & !is.na(tms_segments)
tms_counts_filt <- tms_counts[, tms_keep]
tms_barcodes_filt <- tms_barcodes[tms_keep, ]

message("  TMS after removing DCT/TAL: ", ncol(tms_counts_filt), " cells retained")
message("  TMS cell types kept:")
print(table(tms_barcodes_filt$segment))

# ── 5. Intersect gene sets ──────────────────────────────────────────────────
message("\n▶ Intersecting gene sets…")
shared_genes <- intersect(rownames(tms_counts_filt), rownames(chen_counts))
message("  TMS genes: ", nrow(tms_counts_filt))
message("  Chen genes: ", nrow(chen_counts))
message("  Shared: ", length(shared_genes))

tms_counts_filt <- tms_counts_filt[shared_genes, ]
chen_counts <- chen_counts[shared_genes, ]

# ── 6. Merge into hybrid matrix ─────────────────────────────────────────────
message("\n▶ Merging into hybrid matrix…")

# Build unified metadata
tms_meta_df <- data.frame(
    barcode = tms_barcodes_filt$index,
    cell_type = tms_barcodes_filt$segment,
    donor_id = tms_barcodes_filt$donor_id,
    source = "TMS",
    stringsAsFactors = FALSE
)

chen_meta_df <- data.frame(
    barcode = colnames(chen_counts),
    cell_type = chen_meta$cell_type,
    donor_id = chen_meta$donor_id,
    source = "Chen",
    stringsAsFactors = FALSE
)

hybrid_meta <- rbind(tms_meta_df, chen_meta_df)

# Merge count matrices
hybrid_counts <- cbind(tms_counts_filt, chen_counts)

message("  Hybrid: ", ncol(hybrid_counts), " cells × ", nrow(hybrid_counts), " genes")
message("  Cell type distribution:")
print(table(hybrid_meta$cell_type))
message("  Source distribution:")
print(table(hybrid_meta$source))
message("  Donor distribution:")
print(table(hybrid_meta$donor_id))

# ── 7. Write output ─────────────────────────────────────────────────────────
message("\n▶ Writing hybrid reference to ", OUT_DIR, "…")

# Write matrix.mtx
writeMM(hybrid_counts, file.path(OUT_DIR, "matrix.mtx"))

# Write barcodes.tsv in same format as TMS
# The key columns deconvolution.R needs: index, free_annotation, donor_id
barcodes_out <- data.frame(
    index            = hybrid_meta$barcode,
    free_annotation  = hybrid_meta$cell_type,
    donor_id         = hybrid_meta$donor_id,
    source           = hybrid_meta$source,
    stringsAsFactors = FALSE
)
write.table(barcodes_out, file.path(OUT_DIR, "barcodes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

# Write features.tsv in same format as TMS
# Need: Ensembl ID (row name), feature_name (symbol)
# Map Ensembl back to symbols for feature_name
feature_symbols <- mapIds(org.Mm.eg.db,
    keys = shared_genes,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
)
feature_symbols[is.na(feature_symbols)] <- shared_genes[is.na(feature_symbols)]

features_out <- data.frame(
    feature_name = feature_symbols,
    row.names = shared_genes,
    stringsAsFactors = FALSE
)
write.table(features_out, file.path(OUT_DIR, "features.tsv"),
    sep = "\t", quote = FALSE, col.names = NA
)

message("\n✓ Hybrid reference built successfully!")
message("  Output: ", OUT_DIR)
message("  Files: matrix.mtx, barcodes.tsv, features.tsv")
message("  Total cells: ", ncol(hybrid_counts))
message("  Total genes: ", nrow(hybrid_counts))
