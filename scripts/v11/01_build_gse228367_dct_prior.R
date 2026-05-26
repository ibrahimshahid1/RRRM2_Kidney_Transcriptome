#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "data/results/run_20260526_v11_dct1_phospho_mediation/dct_prior"
qc_dir <- if (length(args) >= 2) args[[2]] else "data/results/run_20260526_v11_dct1_phospho_mediation/external_qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

dct1_path <- "data/external/dct_reference/GSE228367/GSE228367_CONTROL_DCT1.rds.gz"
dct2_path <- "data/external/dct_reference/GSE228367/GSE228367_CONTROL_DCT2.rds.gz"

load_double_gz_rds <- function(path) {
  con <- gzcon(gzfile(path, "rb"))
  on.exit(close(con), add = TRUE)
  readRDS(con)
}

get_rna_data <- function(obj) {
  assays <- attr(obj, "assays")
  if (!"RNA" %in% names(assays)) {
    stop("RNA assay not found")
  }
  attr(assays$RNA, "data")
}

sparse_gene_stats <- function(mat, metadata, subtype) {
  a <- attributes(mat)
  dims <- a$Dim
  genes <- a$Dimnames[[1]]
  cells <- a$Dimnames[[2]]
  i <- a$i + 1L
  p <- a$p
  x <- a$x
  n_genes <- dims[[1]]
  n_cells <- dims[[2]]
  if (length(cells) != nrow(metadata)) {
    stop("Metadata rows do not match sparse matrix columns")
  }

  weighted_tabulate <- function(index, weights, nbins) {
    out <- numeric(nbins)
    if (length(index) == 0) {
      return(out)
    }
    summed <- rowsum(weights, group = index, reorder = FALSE)
    out[as.integer(rownames(summed))] <- as.numeric(summed[, 1])
    out
  }

  gene_sum <- weighted_tabulate(i, x, n_genes)
  gene_detect <- tabulate(i, nbins = n_genes)
  overall <- data.frame(
    gene_symbol = genes,
    subtype = subtype,
    n_cells = n_cells,
    mean_expr = gene_sum / n_cells,
    pct_detected = gene_detect / n_cells,
    stringsAsFactors = FALSE
  )

  col_nnz <- diff(p)
  col_index <- rep.int(seq_len(n_cells), col_nnz)
  reps <- as.character(metadata$Rep)
  rep_levels <- sort(unique(reps))
  rep_tables <- vector("list", length(rep_levels))
  names(rep_tables) <- rep_levels

  for (rep_name in rep_levels) {
    keep_nz <- reps[col_index] == rep_name
    keep_cols <- reps == rep_name
    n_rep_cells <- sum(keep_cols)
    rep_sum <- weighted_tabulate(i[keep_nz], x[keep_nz], n_genes)
    rep_detect <- tabulate(i[keep_nz], nbins = n_genes)
    rep_tables[[rep_name]] <- data.frame(
      gene_symbol = genes,
      subtype = subtype,
      rep = rep_name,
      n_cells = n_rep_cells,
      mean_expr = rep_sum / n_rep_cells,
      pct_detected = rep_detect / n_rep_cells,
      stringsAsFactors = FALSE
    )
  }

  list(overall = overall, by_rep = do.call(rbind, rep_tables))
}

object_inventory <- function(obj, path, label) {
  md <- attr(obj, "meta.data")
  mat <- get_rna_data(obj)
  a <- attributes(mat)
  data.frame(
    dataset = "GSE228367",
    label = label,
    file = path,
    object_class = paste(class(obj), collapse = ";"),
    n_cells = nrow(md),
    n_metadata_columns = ncol(md),
    assay_names = paste(names(attr(obj, "assays")), collapse = ";"),
    rna_n_genes = a$Dim[[1]],
    rna_n_cells = a$Dim[[2]],
    reps = paste(names(table(md$Rep)), table(md$Rep), sep = ":", collapse = ";"),
    diets = paste(names(table(md$Diet)), table(md$Diet), sep = ":", collapse = ";"),
    subclass_l1 = paste(names(table(md$subclass.l1)), table(md$subclass.l1), sep = ":", collapse = ";"),
    class_all = paste(names(table(md$class_all)), table(md$class_all), sep = ":", collapse = ";"),
    stringsAsFactors = FALSE
  )
}

marker_qc <- function(stats, markers) {
  rows <- stats[stats$gene_symbol %in% markers, ]
  rows[order(match(rows$gene_symbol, markers), rows$subtype), ]
}

cat("Loading GSE228367 DCT1...\n")
dct1 <- load_double_gz_rds(dct1_path)
cat("Loading GSE228367 DCT2...\n")
dct2 <- load_double_gz_rds(dct2_path)

inv <- rbind(
  object_inventory(dct1, dct1_path, "CONTROL_DCT1"),
  object_inventory(dct2, dct2_path, "CONTROL_DCT2")
)
write.table(inv, file.path(qc_dir, "gse228367_object_inventory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Extracting sparse gene stats...\n")
dct1_stats <- sparse_gene_stats(get_rna_data(dct1), attr(dct1, "meta.data"), "DCT1")
dct2_stats <- sparse_gene_stats(get_rna_data(dct2), attr(dct2, "meta.data"), "DCT2")

overall <- rbind(dct1_stats$overall, dct2_stats$overall)
by_rep <- rbind(dct1_stats$by_rep, dct2_stats$by_rep)
write.table(overall, file.path(out_dir, "gse228367_gene_stats_overall.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(by_rep, file.path(out_dir, "gse228367_gene_stats_by_rep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

markers <- unique(c(
  "Slc12a3", "Pvalb", "Trpm6", "Klhl3", "Scnn1a", "Trpv5", "Calb1",
  "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Cul3", "Nedd4l", "Sgk1",
  "Kcnj10", "Kcnj16"
))
write.table(marker_qc(overall, markers), file.path(qc_dir, "gse228367_marker_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

wide_overall <- merge(
  dct1_stats$overall[, c("gene_symbol", "mean_expr", "pct_detected", "n_cells")],
  dct2_stats$overall[, c("gene_symbol", "mean_expr", "pct_detected", "n_cells")],
  by = "gene_symbol",
  all = TRUE,
  suffixes = c("_dct1", "_dct2")
)
wide_overall[is.na(wide_overall)] <- 0

rep_wide <- reshape(
  by_rep[, c("gene_symbol", "subtype", "rep", "mean_expr", "pct_detected")],
  idvar = c("gene_symbol", "subtype"),
  timevar = "rep",
  direction = "wide"
)

rep_dct1 <- rep_wide[rep_wide$subtype == "DCT1", ]
rep_dct2 <- rep_wide[rep_wide$subtype == "DCT2", ]
rep_cols <- grep("^mean_expr\\.", names(rep_dct1), value = TRUE)

rep_de <- merge(
  rep_dct1[, c("gene_symbol", rep_cols)],
  rep_dct2[, c("gene_symbol", rep_cols)],
  by = "gene_symbol",
  all = TRUE,
  suffixes = c("_dct1", "_dct2")
)
rep_de[is.na(rep_de)] <- 0

calc_welch <- function(v1, v2) {
  n1 <- length(v1)
  n2 <- length(v2)
  m1 <- mean(v1)
  m2 <- mean(v2)
  s1 <- stats::var(v1)
  s2 <- stats::var(v2)
  se <- sqrt(s1 / n1 + s2 / n2)
  if (!is.finite(se) || se == 0) {
    return(c(delta = m1 - m2, se = se, t = NA, df = NA, p = ifelse(m1 == m2, 1, NA)))
  }
  tval <- (m1 - m2) / se
  df <- (s1 / n1 + s2 / n2)^2 / ((s1 / n1)^2 / (n1 - 1) + (s2 / n2)^2 / (n2 - 1))
  pval <- 2 * stats::pt(-abs(tval), df = df)
  c(delta = m1 - m2, se = se, t = tval, df = df, p = pval)
}

dct1_rep_cols <- paste0(rep_cols, "_dct1")
dct2_rep_cols <- paste0(rep_cols, "_dct2")
welch <- t(vapply(seq_len(nrow(rep_de)), function(row_i) {
  calc_welch(as.numeric(rep_de[row_i, dct1_rep_cols]), as.numeric(rep_de[row_i, dct2_rep_cols]))
}, numeric(5)))
welch <- as.data.frame(welch)

de <- merge(wide_overall, data.frame(gene_symbol = rep_de$gene_symbol, welch, stringsAsFactors = FALSE), by = "gene_symbol", all.x = TRUE)
eps <- 0.01
de$log2_mean_ratio_dct1_vs_dct2 <- log2((de$mean_expr_dct1 + eps) / (de$mean_expr_dct2 + eps))
de$dct1_enrichment_score <- de$mean_expr_dct1 - de$mean_expr_dct2
de$p_value <- de$p
de$fdr <- p.adjust(ifelse(is.na(de$p_value), 1, de$p_value), method = "BH")
de$dct_expression_class <- "unresolved"
de$dct_expression_class[
  de$log2_mean_ratio_dct1_vs_dct2 > 1 &
    de$fdr < 0.05 &
    de$pct_detected_dct1 >= 0.05 &
    de$mean_expr_dct1 >= 0.05
] <- "DCT1_core"
de$dct_expression_class[
  de$log2_mean_ratio_dct1_vs_dct2 < -1 &
    de$fdr < 0.05 &
    de$pct_detected_dct2 >= 0.05 &
    de$mean_expr_dct2 >= 0.05
] <- "DCT2_core"
de$dct_expression_class[
  de$dct_expression_class == "unresolved" &
    de$pct_detected_dct1 >= 0.05 &
    de$pct_detected_dct2 >= 0.05 &
    abs(de$log2_mean_ratio_dct1_vs_dct2) <= 1
] <- "DCT_shared"

de <- de[order(-de$dct1_enrichment_score), ]
write.table(de, file.path(out_dir, "gse228367_dct1_vs_dct2_de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(de[de$dct_expression_class == "DCT1_core", ], file.path(out_dir, "gse228367_dct1_core_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(de[de$dct_expression_class == "DCT2_core", ], file.path(out_dir, "gse228367_dct2_core_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(de[de$dct_expression_class == "DCT_shared", ], file.path(out_dir, "gse228367_dct_shared_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_lines <- c(
  "{",
  sprintf('  "n_genes": %d,', nrow(de)),
  sprintf('  "n_dct1_core": %d,', sum(de$dct_expression_class == "DCT1_core")),
  sprintf('  "n_dct2_core": %d,', sum(de$dct_expression_class == "DCT2_core")),
  sprintf('  "n_dct_shared": %d,', sum(de$dct_expression_class == "DCT_shared")),
  sprintf('  "n_unresolved": %d,', sum(de$dct_expression_class == "unresolved")),
  '  "primary_score": "dct1_enrichment_score = mean log-normalized RNA expression in DCT1 minus DCT2",',
  '  "core_rule": "abs(log2((mean+0.01)/(mean+0.01))) > 1, FDR < 0.05, pct detected >= 0.05, mean expression >= 0.05",',
  '  "statistical_caveat": "Welch tests use three NK replicate-level mean-expression profiles per subtype; this is a reference-prior analysis, not spaceflight evidence."',
  "}"
)
writeLines(summary_lines, file.path(out_dir, "gse228367_dct_prior_summary.json"))
cat("Done. Outputs in", out_dir, "and", qc_dir, "\n")
