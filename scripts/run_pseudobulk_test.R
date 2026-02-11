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

# USER PATHS
bulk_counts_csv <- "data/processed/aligned_outputs/rsem_rRNArm_raw_counts.csv"
bulk_meta_csv <- "data/processed/aligned_outputs/metadata_aligned.tsv"
sc_h5ad_path <- "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"
outdir <- "data/processed/deconvolution/pseudobulk_test"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- HELPER FUNCTIONS ---
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)
detect_id_space <- function(ids) {
    ids <- ids[!is.na(ids)]
    if (length(ids) == 0) {
        return("unknown")
    }
    frac_ens <- mean(grepl("^ENSMUSG", ids))
    if (frac_ens > 0.5) "ensembl" else "symbol"
}
symbols_to_ensembl <- function(symbols) {
    map <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENSEMBL")
    unique(na.omit(map$ENSEMBL))
}
ensembl_to_symbols <- function(ensembl_ids) {
    map <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(ensembl_ids), keytype = "ENSEMBL", columns = "SYMBOL")
    unique(na.omit(map$SYMBOL))
}
resolve_markers_to_matrix <- function(markers, matrix_rownames) {
    space <- detect_id_space(matrix_rownames)
    looks_ens <- grepl("^ENSMUSG", markers)
    if (space == "ensembl") {
        ens <- if (all(looks_ens)) unique(markers) else symbols_to_ensembl(markers)
        return(intersect(ens, matrix_rownames))
    } else if (space == "symbol") {
        sym <- if (all(looks_ens)) ensembl_to_symbols(markers) else unique(markers)
        return(intersect(sym, matrix_rownames))
    } else {
        # warning("Could not detect ID space of matrix rownames.")
        return(character(0))
    }
}
dedup_sce_genes_sum <- function(sce_obj) {
    g <- sub("\\.\\d+$", "", rownames(sce_obj))
    if (anyDuplicated(g) == 0) {
        rownames(sce_obj) <- g
        return(sce_obj)
    }
    message("SCE has duplicate gene IDs after stripping versions. Collapsing by sum...")
    f <- factor(g, levels = unique(g))
    grp <- as.integer(f)
    G <- sparseMatrix(i = grp, j = seq_along(grp), x = 1, dims = c(nlevels(f), length(grp)))
    new_counts <- G %*% counts(sce_obj)
    out <- SingleCellExperiment(assays = list(counts = new_counts), colData = colData(sce_obj))
    rownames(out) <- levels(f)
    out
}
pick_within_group_markers <- function(sce_obj, types, clusters_col = "segment_use", top_n = 120, min_mean = 0.5, lfc_min = 1.0) {
    X <- assay(sce_obj, "counts")
    cl <- as.character(colData(sce_obj)[[clusters_col]])
    keep <- cl %in% types
    if (sum(keep) == 0) {
        return(character(0))
    }
    Xg <- X[, keep, drop = FALSE]
    clg <- factor(cl[keep])
    lib <- Matrix::colSums(Xg)
    Xg_cpm <- t(t(Xg) / (lib + 1e-12)) * 1e6
    Xg_log <- log1p(Xg_cpm)
    means <- sapply(levels(clg), function(t) Matrix::rowMeans(Xg_log[, clg == t, drop = FALSE]))
    rownames(means) <- rownames(Xg)
    markers <- character(0)
    for (t in colnames(means)) {
        mu_t <- means[, t]
        mu_o <- rowMeans(means[, setdiff(colnames(means), t), drop = FALSE])
        lfc <- mu_t - mu_o
        ok <- is.finite(lfc) & (mu_t > min_mean) & (lfc > lfc_min)
        if (sum(ok) == 0) next
        top <- names(sort(lfc[ok], decreasing = TRUE))[seq_len(min(top_n, sum(ok)))]
        markers <- unique(c(markers, top))
    }
    markers
}
pick_type_markers <- function(sce_obj, target, types, clusters_col = "segment_use", top_n = 50, min_mean = 0.1, lfc_min = 0.5) {
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
normalize_rows_to_one <- function(M, eps = 1e-12) {
    rs <- rowSums(M)
    rs[rs < eps] <- 1
    M / rs
}
rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
rdirichlet_base <- function(n, alpha) {
    k <- length(alpha)
    X <- matrix(0, nrow = n, ncol = k)
    for (j in seq_len(k)) X[, j] <- rgamma(n, shape = alpha[j], rate = 1)
    X <- X / rowSums(X)
    X
}
make_pseudobulk_exact <- function(sce_obj, segment_col = "segment_use", donor_col = "donor_id", n_mixtures = 30, cells_per_mixture = 2000, depth_scale = 1.0, target_depth = NULL, use_multinomial_depth = TRUE, dct_range = c(0.00, 0.25), pt_range = c(0.40, 0.75), tal_range = c(0.05, 0.25), cd_range = c(0.00, 0.15), other_tubule_range = c(0.00, 0.10), tubule_total_range = c(0.50, 0.85), glom_total_range = c(0.05, 0.20), immune_frac_range = c(0.01, 0.10), stroma_frac_range = c(0.02, 0.10), seed = 1) {
    set.seed(seed)
    seg <- factor(as.character(colData(sce_obj)[[segment_col]]))
    donors <- unique(as.character(colData(sce_obj)[[donor_col]]))
    donors <- donors[!is.na(donors)]
    seg_levels <- levels(seg)
    genes <- rownames(sce_obj)
    n_genes <- length(genes)

    glom_segs <- intersect(c("Podocyte", "Endothelial", "Mesangial"), seg_levels)
    tubule_segs <- intersect(c("PT", "TAL_LOH", "DCT", "CD", "Other"), seg_levels)
    immune_segs <- intersect(c("Immune"), seg_levels)
    stroma_segs <- intersect(c("Fibroblast"), seg_levels)

    cell_idx <- seq_len(ncol(sce_obj))
    idx_by_donor_seg <- list()
    for (d in donors) idx_by_donor_seg[[d]] <- split(cell_idx[colData(sce_obj)[[donor_col]] == d], seg[colData(sce_obj)[[donor_col]] == d])

    # Ensure robust donor list
    if (length(donors) == 0) stop("No valid donors found!")

    pb_mat <- matrix(0, nrow = n_genes, ncol = n_mixtures, dimnames = list(genes, paste0("pb", seq_len(n_mixtures))))
    P_rna_true <- matrix(0, nrow = n_mixtures, ncol = length(seg_levels), dimnames = list(paste0("pb", seq_len(n_mixtures)), seg_levels))

    for (i in seq_len(n_mixtures)) {
        d <- sample(donors, 1)

        tubule_total <- runif(1, tubule_total_range[1], tubule_total_range[2])
        glom_total <- runif(1, glom_total_range[1], glom_total_range[2])
        immune_total <- runif(1, immune_frac_range[1], immune_frac_range[2])
        stroma_total <- runif(1, stroma_frac_range[1], stroma_frac_range[2])

        coarse_raw <- c(Tubule = tubule_total, Glomerular = glom_total, Immune = immune_total, Stroma = stroma_total)
        coarse_norm <- coarse_raw / sum(coarse_raw)

        tubule_within <- c(DCT = runif(1, dct_range[1], dct_range[2]), PT = runif(1, pt_range[1], pt_range[2]), TAL_LOH = runif(1, tal_range[1], tal_range[2]), CD = runif(1, cd_range[1], cd_range[2]), Other = runif(1, other_tubule_range[1], other_tubule_range[2]))
        tubule_within <- tubule_within[names(tubule_within) %in% tubule_segs]
        if (sum(tubule_within) > 0) tubule_within <- tubule_within / sum(tubule_within)

        glom_within <- if (length(glom_segs) > 0) rdirichlet_base(1, rep(2, length(glom_segs)))[1, ] else numeric(0)
        names(glom_within) <- glom_segs

        target_props <- setNames(rep(0, length(seg_levels)), seg_levels)
        for (s in names(tubule_within)) if (s %in% seg_levels) target_props[s] <- coarse_norm["Tubule"] * tubule_within[s]
        for (s in names(glom_within)) if (s %in% seg_levels) target_props[s] <- coarse_norm["Glomerular"] * glom_within[s]
        for (s in immune_segs) if (s %in% seg_levels) target_props[s] <- coarse_norm["Immune"] / length(immune_segs)
        for (s in stroma_segs) if (s %in% seg_levels) target_props[s] <- coarse_norm["Stroma"] / length(stroma_segs)

        if (sum(target_props) > 0) target_props <- target_props / sum(target_props)

        n_per_seg <- as.integer(rmultinom(1, size = cells_per_mixture, prob = target_props))
        names(n_per_seg) <- seg_levels

        umi_by_seg <- setNames(rep(0, length(seg_levels)), seg_levels)
        x_total <- numeric(n_genes)
        names(x_total) <- genes

        for (s in seg_levels) {
            if (n_per_seg[s] <= 0) next
            pool <- idx_by_donor_seg[[d]][[s]]
            if (is.null(pool) || length(pool) == 0) next # Donor missing segment

            chosen_s <- sample(pool, size = n_per_seg[s], replace = TRUE)
            xs <- Matrix::rowSums(counts(sce_obj)[, chosen_s, drop = FALSE])
            umi_by_seg[s] <- sum(xs)
            x_total <- x_total + as.numeric(xs)
        }

        if (sum(umi_by_seg) > 0) P_rna_true[i, ] <- umi_by_seg / sum(umi_by_seg)
        pb_mat[, i] <- x_total
    }

    list(pb_mat = pb_mat, P_rna_true = P_rna_true)
}

run_music_on_pseudobulk <- function(pb_mat, sce_obj, group.markers, clusters.type, normalize_flag = FALSE) {
    storage.mode(pb_mat) <- "numeric"
    common <- intersect(rownames(pb_mat), rownames(sce_obj))
    pb_mat2 <- pb_mat[common, , drop = FALSE]
    sce2b <- sce_obj[common, ]

    MuSiC::music_prop.cluster(
        bulk.mtx = pb_mat2, sc.sce = sce2b, group.markers = group.markers,
        groups = "clusterType", clusters = "segment_use", samples = "donor_id",
        clusters.type = clusters.type, normalize = normalize_flag, centered = FALSE, verbose = FALSE
    )
}

evaluate_pseudobulk <- function(P_true, P_hat) {
    common_segs <- intersect(colnames(P_true), colnames(P_hat))
    T <- P_true[, common_segs, drop = FALSE]
    H <- P_hat[, common_segs, drop = FALSE]
    stats <- data.frame(segment = common_segs, spearman_rho = NA_real_, rmse = NA_real_)
    for (j in seq_along(common_segs)) {
        s <- common_segs[j]
        stats$spearman_rho[j] <- suppressWarnings(cor(T[, s], H[, s], method = "spearman"))
        stats$rmse[j] <- rmse(T[, s], H[, s])
    }
    list(stats = stats)
}

# --- MAIN EXECUTION ---
message("Loading single-cell reference from Raw MTX...")
raw_dir <- "data/processed/deconvolution/raw_ref"
mtx_path <- file.path(raw_dir, "matrix.mtx")
obs_path <- file.path(raw_dir, "barcodes.tsv")
var_path <- file.path(raw_dir, "features.tsv")

raw_mat <- readMM(mtx_path)
obs <- data.frame(fread(obs_path, header = TRUE), row.names = 1, check.names = FALSE)
var <- data.frame(fread(var_path, header = TRUE), row.names = 1, check.names = FALSE)

# Transpose check
if (ncol(raw_mat) == nrow(var)) raw_mat <- t(raw_mat)
colnames(raw_mat) <- rownames(obs)
rownames(raw_mat) <- rownames(var)

sce <- SingleCellExperiment(assays = list(counts = raw_mat), colData = obs)
sce <- dedup_sce_genes_sum(sce)

# Label processing
label_cols <- intersect(c("free_annotation", "cell_type", "subtissue"), colnames(colData(sce)))
score_col <- function(col) sum(grepl("podo|glom", as.character(colData(sce)[[col]]), ignore.case = TRUE))
scores <- sapply(label_cols, score_col)
best <- if ("free_annotation" %in% names(scores)) "free_annotation" else names(which.max(scores))
colData(sce)$cell_type_use <- factor(as.character(colData(sce)[[best]]))

ct <- as.character(colData(sce)$cell_type_use)
ct[tolower(trimws(ct)) %in% c("", "nan", "na", "null", "unknown")] <- NA
sce <- sce[, !is.na(ct)]
colData(sce)$cell_type_use <- factor(ct[!is.na(ct)])

# Mapping
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
    "Other"
}
seg <- vapply(as.character(colData(sce)$cell_type_use), segment_from_ct, character(1))
colData(sce)$segment_use <- factor(seg)
sce <- sce[, !is.na(colData(sce)$segment_use)]
colData(sce)$segment_use <- droplevels(colData(sce)$segment_use)
colData(sce)$segment_use <- as.character(colData(sce)$segment_use)

# EXP_REMOVE_OTHER
sce2 <- sce[, colData(sce)$segment_use != "Other"]
colData(sce2)$segment_use <- droplevels(factor(colData(sce2)$segment_use))

seg <- as.character(colData(sce2)$segment_use)
colData(sce2)$clusterType <- factor(ifelse(seg %in% c("Podocyte", "Endothelial", "Mesangial"), "Glomerular",
    ifelse(seg == "Immune", "Immune",
        ifelse(seg == "Fibroblast", "Stroma",
            ifelse(seg == "PT", "Proximal", "Distal")
        )
    )
))

clusters.type <- list(
    Glomerular = c("Podocyte", "Endothelial", "Mesangial"),
    Proximal = c("PT"),
    Distal = c("TAL_LOH", "DCT", "CD"),
    Immune = c("Immune"),
    Stroma = c("Fibroblast")
)

# Marker Selection
message("Selecting markers...")
tub_types <- intersect(c("PT", "TAL_LOH", "DCT", "CD"), unique(colData(sce2)$segment_use))
mk_pt <- pick_type_markers(sce2, "PT", tub_types, top_n = 50)
mk_tal_final <- unique(c(resolve_markers_to_matrix(c("Umod", "Slc12a1", "Cldn10", "Cldn16", "Kcnj1", "Bsnd"), rownames(sce2)), pick_type_markers(sce2, "TAL_LOH", tub_types, top_n = 15)))
mk_dct_final <- unique(c(resolve_markers_to_matrix(c("Slc12a3", "Pvalb", "Calb1", "Trpv5", "Kl", "Wnk4"), rownames(sce2)), pick_type_markers(sce2, "DCT", tub_types, top_n = 15)))
mk_cd_final <- unique(c(resolve_markers_to_matrix(c("Aqp2", "Aqp3", "Fxyd4", "Hsd11b2", "Scnn1g", "Krt8", "Krt18"), rownames(sce2)), pick_type_markers(sce2, "CD", tub_types, top_n = 15)))

group.markers <- list(
    Glomerular = pick_within_group_markers(sce2, clusters.type$Glomerular, top_n = 50),
    Proximal = unique(c(resolve_markers_to_matrix(c("Slc34a1", "Lrp2", "Aqp1", "Slc22a6", "Slc5a2"), rownames(sce2)), mk_pt)),
    Distal = unique(c(resolve_markers_to_matrix(c("Umod", "Slc12a1", "Slc12a3", "Pvalb", "Calb1", "Aqp2"), rownames(sce2)), mk_tal_final, mk_dct_final, mk_cd_final)),
    Immune = NULL,
    Stroma = NULL
)

message("Running Pseudobulk Test with 20 mixtures...")
pb <- make_pseudobulk_exact(sce2, n_mixtures = 20, cells_per_mixture = 2000)
if ("P_rna_true" %in% names(pb)) pb$P_rna <- pb$P_rna_true

res_T <- run_music_on_pseudobulk(pb$pb_mat, sce2, group.markers, clusters.type, normalize_flag = TRUE)
eval_T <- evaluate_pseudobulk(pb$P_rna[rownames(res_T$Est.prop.weighted.cluster), ], res_T$Est.prop.weighted.cluster)

print(eval_T$stats)
