#!/usr/bin/env Rscript
# ---------------------------
# DCT Marker Diagnostic Script
# ---------------------------
# Run this AFTER run_deconvolution.R (requires sce2, bulk_mat2 to be loaded)

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

message("\n", paste(rep("=", 70), collapse = ""))
message("DCT MARKER DIAGNOSTIC")
message(paste(rep("=", 70), collapse = ""))

# Classic DCT markers from literature
dct_markers_symbols <- c(
  "Slc12a3",   # NCC sodium-chloride cotransporter (THE canonical DCT marker)
  "Pvalb",     # Parvalbumin (DCT1 marker)
  "Trpv5",     # Calcium channel (DCT2/CNT)
  "Calb1",     # Calbindin (DCT2)
  "Fxyd2",     # FXYD domain containing 2
  "Trpm6",     # Magnesium channel
  "Klk1",      # Kallikrein 1 (DCT)
  "Slc8a1",    # NCX1 sodium-calcium exchanger
  "Wnk1",      # With-no-lysine kinase 1
  "Pth1r"      # PTH receptor (PT/DCT)
)

# ---------------------------
# 1. Check marker presence in reference and bulk
# ---------------------------
message("\n--- 1. Marker Presence Check ---")

# Check which ID space we're in
ref_space <- if (mean(grepl("^ENSMUSG", rownames(sce2))) > 0.5) "Ensembl" else "Symbol"
bulk_space <- if (mean(grepl("^ENSMUSG", rownames(bulk_mat2))) > 0.5) "Ensembl" else "Symbol"

message(sprintf("Reference ID space: %s", ref_space))
message(sprintf("Bulk ID space: %s", bulk_space))

# If Ensembl, try to find markers anyway (maybe they're symbols in var)
in_ref <- intersect(dct_markers_symbols, rownames(sce2))
in_bulk <- intersect(dct_markers_symbols, rownames(bulk_mat2))

message(sprintf("\nDCT markers found in reference (%d/%d): %s", 
                length(in_ref), length(dct_markers_symbols),
                if (length(in_ref) > 0) paste(in_ref, collapse = ", ") else "NONE"))
message(sprintf("DCT markers found in bulk (%d/%d): %s",
                length(in_bulk), length(dct_markers_symbols),
                if (length(in_bulk) > 0) paste(in_bulk, collapse = ", ") else "NONE"))

# If Ensembl, try mapping
if (ref_space == "Ensembl" && length(in_ref) == 0) {
  message("\nAttempting Ensembl mapping for reference...")
  if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    
    map <- AnnotationDbi::select(org.Mm.eg.db,
                                  keys = dct_markers_symbols,
                                  keytype = "SYMBOL",
                                  columns = "ENSEMBL")
    ens_ids <- na.omit(map$ENSEMBL)
    in_ref <- intersect(ens_ids, rownames(sce2))
    message(sprintf("After Ensembl mapping: %d/%d markers found in reference",
                    length(in_ref), length(dct_markers_symbols)))
    
    # Create mapping for display
    marker_map <- setNames(map$ENSEMBL[!is.na(map$ENSEMBL)], 
                           map$SYMBOL[!is.na(map$ENSEMBL)])
  }
}

if (bulk_space == "Ensembl" && length(in_bulk) == 0) {
  message("\nAttempting Ensembl mapping for bulk...")
  if (exists("marker_map")) {
    ens_ids <- unname(marker_map)
    in_bulk <- intersect(ens_ids, rownames(bulk_mat2))
    message(sprintf("After Ensembl mapping: %d/%d markers found in bulk",
                    length(in_bulk), length(dct_markers_symbols)))
  }
}

# ---------------------------
# 2. Expression specificity in reference
# ---------------------------
if (length(in_ref) > 0) {
  message("\n--- 2. DCT Marker Specificity in Reference ---")
  
  X <- assay(sce2, "counts")
  seg <- as.character(colData(sce2)$segment_use)
  
  # Mean expression per segment
  segments <- unique(seg)
  expr_by_seg <- sapply(segments, function(s) {
    cells <- seg == s
    Matrix::rowMeans(X[in_ref, cells, drop = FALSE])
  })
  rownames(expr_by_seg) <- in_ref
  
  # For display, try to show symbol names
  if (exists("marker_map")) {
    rev_map <- setNames(names(marker_map), marker_map)
    display_names <- sapply(in_ref, function(x) {
      if (x %in% names(rev_map)) paste0(rev_map[x], " (", x, ")") else x
    })
    rownames(expr_by_seg) <- display_names
  }
  
  message("\nMean expression by segment:")
  print(round(expr_by_seg, 2))
  
  # Calculate specificity (DCT / max other)
  message("\n--- DCT Marker Specificity Ratios ---")
  message("(DCT expression / max other segment expression)")
  
  for (i in seq_len(nrow(expr_by_seg))) {
    gene <- rownames(expr_by_seg)[i]
    dct_expr <- expr_by_seg[i, "DCT"]
    other_expr <- expr_by_seg[i, setdiff(colnames(expr_by_seg), "DCT")]
    max_other <- max(other_expr)
    which_max <- names(which.max(other_expr))
    ratio <- dct_expr / (max_other + 0.01)
    
    status <- if (ratio > 2) "✓ GOOD" else if (ratio > 1) "~ OK" else "✗ POOR"
    message(sprintf("  %s: DCT=%.2f, MaxOther=%.2f (%s), Ratio=%.1fx %s",
                    gene, dct_expr, max_other, which_max, ratio, status))
  }
  
  # ---------------------------
  # 3. Bulk expression of DCT markers
  # ---------------------------
  message("\n--- 3. Bulk Expression of DCT Markers ---")
  
  common_markers <- intersect(in_ref, in_bulk)
  if (length(common_markers) > 0) {
    bulk_expr <- bulk_mat2[common_markers, , drop = FALSE]
    
    message(sprintf("\nBulk mean counts for DCT markers (n=%d):", length(common_markers)))
    print(round(rowMeans(bulk_expr), 1))
    
    # Compute bulk DCT score
    bulk_dct_score <- colMeans(log1p(bulk_expr))
    
    message("\n--- 4. DCT Bulk Score vs Estimated Proportion ---")
    if ("DCT" %in% colnames(prop_seg)) {
      ct <- cor.test(bulk_dct_score, prop_seg[, "DCT"], method = "spearman")
      message(sprintf("Spearman correlation: %.3f (p=%.3e)", ct$estimate, ct$p.value))
      
      if (ct$estimate > 0.3) {
        message("✓ Good correlation - literature markers validate DCT estimates")
      } else if (ct$estimate > 0.1) {
        message("~ Weak correlation - some signal but noisy")
      } else {
        message("✗ Poor correlation - DCT estimates may not be reliable")
      }
    }
  } else {
    message("No common DCT markers between reference and bulk!")
  }
  
} else {
  message("\n✗ CRITICAL: No DCT markers found in reference!")
  message("This strongly suggests the reference is inadequate for DCT deconvolution.")
}

# ---------------------------
# 4. What markers IS MuSiC using for DCT?
# ---------------------------
message("\n--- 5. What Markers is MuSiC Using for DCT? ---")

if (exists("segment_markers") && "DCT" %in% names(segment_markers)) {
  dct_music_markers <- segment_markers$DCT
  message(sprintf("MuSiC derived %d DCT markers:", length(dct_music_markers)))
  
  # Show first 10
  message("Top 10 derived markers:")
  print(head(dct_music_markers, 10))
  
  # Check overlap with literature markers
  overlap <- intersect(dct_music_markers, in_ref)
  message(sprintf("\nOverlap with literature DCT markers: %d/%d", 
                  length(overlap), length(in_ref)))
  if (length(overlap) > 0) {
    message("Overlapping: ", paste(overlap, collapse = ", "))
  }
} else {
  message("segment_markers not found - run validation section first")
}

# ---------------------------
# 5. Cell count check
# ---------------------------
message("\n--- 6. Cell Count Check ---")
seg_counts <- table(colData(sce2)$segment_use)
message("\nCells per segment:")
print(sort(seg_counts, decreasing = TRUE))

dct_count <- seg_counts["DCT"]
total_cells <- sum(seg_counts)
message(sprintf("\nDCT: %d cells (%.1f%% of reference)", dct_count, 100*dct_count/total_cells))

if (dct_count < 500) {
  message("⚠ WARNING: Low DCT cell count may limit deconvolution accuracy")
}
if (dct_count < 200) {
  message("✗ CRITICAL: Very low DCT cell count - consider using a different atlas")
}

message("\n", paste(rep("=", 70), collapse = ""))
message("DIAGNOSTIC COMPLETE")
message(paste(rep("=", 70), collapse = ""))
