#!/usr/bin/env Rscript
# src/networks/wgcna_analysis.R

suppressPackageStartupMessages({
  library(WGCNA)
  library(optparse)
})

# Disable WGCNA threading to avoid issues in pipeline context
disableWGCNAThreads()

options(stringsAsFactors = FALSE)

# CLI args

option_list <- list(
  make_option("--rtech", type = "character",
              help = "Path to Rtech.tsv.gz (genes × samples, tech-corrected)"),
  make_option("--meta", type = "character",
              help = "Path to meta_phase1.tsv.gz"),
  make_option("--id_map", type = "character", default = "",
              help = "Path to id_map.tsv (Ensembl → symbol)"),
  make_option("--gene_sets", type = "character", default = "",
              help = "Path to gene_sets.yaml for module enrichment"),
  make_option("--outdir", type = "character",
              help = "Output directory"),
  make_option("--max_genes", type = "integer", default = 5000,
              help = "Top genes by variance [default %default]"),
  make_option("--power", type = "integer", default = 0,
              help = "Soft threshold power (0 = auto-select) [default %default]"),
  make_option("--min_module_size", type = "integer", default = 30,
              help = "Minimum module size [default %default]"),
  make_option("--merge_cut", type = "double", default = 0.25,
              help = "Module merge cut height [default %default]"),
  make_option("--n_pres_perms", type = "integer", default = 200,
              help = "Permutations for module preservation [default %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed [default %default]")
)

args <- parse_args(OptionParser(option_list = option_list))

set.seed(args$seed)
outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("WGCNA Module Analysis for RRRM-2\n")
cat("============================================================\n\n")

# Load data

cat("Loading expression:", args$rtech, "\n")
rtech <- read.csv(gzfile(args$rtech), sep = "\t", row.names = 1, check.names = FALSE)
cat("  Dimensions:", nrow(rtech), "genes ×", ncol(rtech), "samples\n")

cat("Loading metadata:", args$meta, "\n")
meta <- read.csv(gzfile(args$meta), sep = "\t", check.names = FALSE)

# Find sample column
sample_col <- NA
for (col in c("Sample Name (raw_counts_colname)", "Sample Name", "sample")) {
  if (col %in% colnames(meta)) { sample_col <- col; break }
}
if (is.na(sample_col)) sample_col <- colnames(meta)[1]
rownames(meta) <- meta[[sample_col]]

# Normalize labels
if ("EnvGroup" %in% colnames(meta)) {
  meta$EnvGroup[meta$EnvGroup == "HGC"] <- "GC"
  meta$EnvGroup[meta$EnvGroup == "VGC"] <- "VIV"
}
if ("Age" %in% colnames(meta)) {
  meta$Age[meta$Age %in% c("Young", "young", "Yng", "YOUNG")] <- "YNG"
  meta$Age[meta$Age %in% c("old")] <- "OLD"
}
if ("Arm" %in% colnames(meta)) {
  meta$Arm[meta$Arm %in% c("ISS", "ISST", "ISS_T")] <- "ISS-T"
}

# Align samples
common <- intersect(colnames(rtech), rownames(meta))
rtech <- rtech[, common]
meta <- meta[common, ]
cat("  Aligned:", length(common), "samples\n")

# Filter to FLT + GC only

keep <- meta$EnvGroup %in% c("FLT", "GC")
cat("\nFiltering to FLT + GC only:", sum(keep), "samples\n")
rtech_fg <- rtech[, keep]
meta_fg <- meta[keep, ]
cat("  FLT:", sum(meta_fg$EnvGroup == "FLT"), ", GC:", sum(meta_fg$EnvGroup == "GC"), "\n")
cat("  Age: ", paste(table(meta_fg$Age), collapse = "/"), "\n")
cat("  Arm: ", paste(table(meta_fg$Arm), collapse = "/"), "\n")

# Gene selection: top N by variance

gene_var <- apply(rtech_fg, 1, var, na.rm = TRUE)
gene_var <- sort(gene_var, decreasing = TRUE)
n_select <- min(args$max_genes, length(gene_var))
top_genes <- names(gene_var)[1:n_select]
cat("\nTop", n_select, "genes by variance selected\n")
cat("  Variance range: [", round(gene_var[n_select], 4), ",",
    round(gene_var[1], 4), "]\n")

# WGCNA expects samples × genes
datExpr <- as.data.frame(t(rtech_fg[top_genes, ]))

# QC
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  cat("  Removing", sum(!gsg$goodGenes), "bad genes,",
      sum(!gsg$goodSamples), "bad samples\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  meta_fg <- meta_fg[gsg$goodSamples, ]
}
cat("  Final matrix:", nrow(datExpr), "samples ×", ncol(datExpr), "genes\n")

# Sample clustering QC

pdf(file.path(outdir, "sample_dendrogram.pdf"), width = 12, height = 6)
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering (FLT + GC)", xlab = "", sub = "",
     labels = paste(meta_fg$EnvGroup, meta_fg$Age, meta_fg$Arm, sep = "_"),
     cex = 0.6)
dev.off()


cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Analysis A: Module Discovery\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

powers <- c(1:20, seq(22, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed",
                         corFnc = "bicor",
                         corOptions = list(maxPOutliers = 0.05),
                         verbose = 0)

# Save soft threshold diagnostics
pdf(file.path(outdir, "soft_threshold_diagnostics.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold power", ylab = "Scale-free topology fit (signed R²)",
     type = "b", pch = 20, main = "Scale-free fit index")
abline(h = 0.80, col = "red", lty = 2)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = sft$fitIndices[, 1], cex = 0.7, pos = 3)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold power", ylab = "Mean connectivity",
     type = "b", pch = 20, main = "Mean connectivity")
dev.off()

write.csv(sft$fitIndices, file.path(outdir, "soft_threshold_fit.csv"), row.names = FALSE)

# Choose power
chosen_power <- args$power
if (chosen_power == 0) {
  # Auto: first power where signed R² > 0.80, or fallback to recommended
  fit_vals <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  above_thresh <- which(fit_vals > 0.80)
  if (length(above_thresh) > 0) {
    chosen_power <- sft$fitIndices[above_thresh[1], 1]
  } else {
    # Conservative fallback: power 12 for signed networks
    chosen_power <- 12
    cat("  WARNING: No power achieved R² > 0.80, using conservative default = 12\n")
  }
}
cat("  Chosen soft threshold power:", chosen_power, "\n")

# Build modules
cat("  Building modules (blockwiseModules)...\n")
net <- blockwiseModules(
  datExpr,
  power = chosen_power,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.05,
  minModuleSize = args$min_module_size,
  mergeCutHeight = args$merge_cut,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,       # skip TOM save to save disk
  verbose = 0
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
MEs <- orderMEs(MEs)
n_modules <- length(unique(moduleColors)) - ("grey" %in% moduleColors)

cat("  Modules found:", n_modules, "(+ grey/unassigned)\n")
cat("  Module sizes:\n")
mod_table <- sort(table(moduleColors), decreasing = TRUE)
print(mod_table)

# Save module assignments
gene_names <- colnames(datExpr)
mod_df <- data.frame(
  gene = gene_names,
  module_num = net$colors,
  module_color = moduleColors,
  stringsAsFactors = FALSE
)

# Add gene symbols if id_map available
if (nchar(args$id_map) > 0 && file.exists(args$id_map)) {
  id_map <- read.csv(args$id_map, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  eid_col <- grep("ensembl", colnames(id_map), ignore.case = TRUE, value = TRUE)[1]
  sym_col <- grep("symbol|mgi", colnames(id_map), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(eid_col) && !is.na(sym_col)) {
    sym_map <- setNames(id_map[[sym_col]], id_map[[eid_col]])
    mod_df$symbol <- sym_map[mod_df$gene]
  }
}

write.csv(mod_df, file.path(outdir, "module_assignments.csv"), row.names = FALSE)
write.csv(as.data.frame(MEs), file.path(outdir, "module_eigengenes.csv"), row.names = TRUE)

# Save per-module gene lists
mod_list_dir <- file.path(outdir, "module_gene_lists")
dir.create(mod_list_dir, showWarnings = FALSE)
for (mc in unique(moduleColors)) {
  genes_in_mod <- gene_names[moduleColors == mc]
  writeLines(genes_in_mod, file.path(mod_list_dir, paste0(mc, ".txt")))
}

# Dendrogram + module colors
pdf(file.path(outdir, "dendrogram_modules.pdf"), width = 14, height = 7)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                   "Module colors", dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05,
                   main = "Gene dendrogram and module colors")
dev.off()



cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Analysis B: Module-Trait Association\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

trait_df <- data.frame(
  Flight = ifelse(meta_fg$EnvGroup == "FLT", 1, 0),
  AgeOld = ifelse(meta_fg$Age == "OLD", 1, 0),
  ArmLAR = ifelse(meta_fg$Arm == "LAR", 1, 0),
  row.names = rownames(meta_fg)
)

# Full factorial model for each module eigengene
trait_results <- list()
for (me_col in colnames(MEs)) {
  df <- data.frame(ME = MEs[, me_col], trait_df)
  fit <- lm(ME ~ Flight * AgeOld * ArmLAR, data = df)
  coefs <- summary(fit)$coefficients
  coefs_df <- as.data.frame(coefs)
  coefs_df$term <- rownames(coefs)
  coefs_df$module <- me_col
  trait_results[[me_col]] <- coefs_df
}

trait_all <- do.call(rbind, trait_results)
rownames(trait_all) <- NULL
colnames(trait_all)[1:4] <- c("estimate", "std_error", "t_value", "p_value")

# BH correction within each term
terms <- unique(trait_all$term)
for (tm in terms) {
  idx <- trait_all$term == tm
  trait_all$q_value[idx] <- p.adjust(trait_all$p_value[idx], method = "BH")
}

write.csv(trait_all, file.path(outdir, "module_trait_association.csv"), row.names = FALSE)

# Summary of significant associations
sig_traits <- trait_all[trait_all$q_value < 0.05 & trait_all$term != "(Intercept)", ]
if (nrow(sig_traits) > 0) {
  cat("  Significant associations (q < 0.05):\n")
  for (i in 1:nrow(sig_traits)) {
    cat(sprintf("    %s ~ %s: est=%.3f, q=%.4f\n",
                sig_traits$module[i], sig_traits$term[i],
                sig_traits$estimate[i], sig_traits$q_value[i]))
  }
} else {
  cat("  No module-trait associations at q < 0.05\n")
}

# Heatmap of module-trait correlations
modTraitCor <- cor(MEs, as.matrix(trait_df), use = "pairwise.complete.obs")
modTraitP <- corPvalueStudent(modTraitCor, nrow(datExpr))

pdf(file.path(outdir, "module_trait_heatmap.pdf"), width = 8, height = 10)
textMatrix <- paste0(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")")
dim(textMatrix) <- dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor,
               xLabels = colnames(trait_df),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()



cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Analysis C: Module Preservation (GC → FLT)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

gc_samples  <- meta_fg$EnvGroup == "GC"
flt_samples <- meta_fg$EnvGroup == "FLT"
cat("  GC samples:", sum(gc_samples), ", FLT samples:", sum(flt_samples), "\n")

multiExpr <- list(
  GC  = list(data = as.data.frame(datExpr[gc_samples, ])),
  FLT = list(data = as.data.frame(datExpr[flt_samples, ]))
)

multiColor <- list(GC = moduleColors)

cat("  Running modulePreservation (", args$n_pres_perms, " permutations)...\n")
# NOTE: Using cor (Pearson) for preservation instead of bicor.
pres <- modulePreservation(
  multiExpr, multiColor,
  referenceNetworks = 1,
  nPermutations = args$n_pres_perms,
  networkType = "signed",
  verbose = 0
)

# Extract preservation statistics
pres_stats <- pres$preservation$Z$ref.GC$inColumnsAlsoPresentIn.FLT
pres_obs <- pres$preservation$observed$ref.GC$inColumnsAlsoPresentIn.FLT

cat("  Available Z columns:", paste(colnames(pres_stats), collapse = ", "), "\n")

# Build summary table with safe column access
get_col <- function(df, patterns) {
  for (p in patterns) {
    matches <- grep(p, colnames(df), value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) return(df[, matches[1]])
  }
  return(rep(NA, nrow(df)))
}

pres_df <- data.frame(
  module_color = rownames(pres_stats),
  moduleSize = get_col(pres_obs, c("moduleSize")),
  Zsummary = get_col(pres_stats, c("Zsummary.pres", "Zsummary")),
  Zdensity = get_col(pres_stats, c("Zdensity.pres", "Zdensity")),
  Zconnectivity = get_col(pres_stats, c("Zconnectivity.pres", "Zconnectivity")),
  medianRank = get_col(pres_stats, c("medianRank.pres", "medianRank")),
  stringsAsFactors = FALSE
)

# Interpretation labels
pres_df$preservation <- ifelse(is.na(pres_df$Zsummary), "unknown",
                        ifelse(pres_df$Zsummary > 10, "strong",
                        ifelse(pres_df$Zsummary > 2, "moderate", "none")))

# Remove gold (WGCNA internal reference module)
pres_df <- pres_df[!pres_df$module_color %in% c("gold"), ]

write.csv(pres_df, file.path(outdir, "module_preservation.csv"), row.names = FALSE)

cat("  Preservation results:\n")
for (i in 1:nrow(pres_df)) {
  cat(sprintf("    %-15s (n=%3d): Zsummary=%6.2f  [%s]\n",
              pres_df$module_color[i],
              ifelse(is.na(pres_df$moduleSize[i]), 0, pres_df$moduleSize[i]),
              ifelse(is.na(pres_df$Zsummary[i]), 0, pres_df$Zsummary[i]),
              pres_df$preservation[i]))
}

# Preservation plot
pdf(file.path(outdir, "module_preservation_plot.pdf"), width = 8, height = 6)
plot_df <- pres_df[pres_df$module_color != "grey", ]
if (nrow(plot_df) > 0) {
  barplot(plot_df$Zsummary,
          names.arg = plot_df$module_color,
          col = plot_df$module_color,
          las = 2, ylab = "Zsummary.pres",
          main = "Module preservation: GC reference → FLT test")
  abline(h = 10, col = "darkgreen", lty = 2)
  abline(h = 2, col = "red", lty = 2)
  legend("topright", c("strong (>10)", "moderate (2-10)", "none (<2)"),
         lty = 2, col = c("darkgreen", "orange", "red"), cex = 0.8)
}
dev.off()



cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Analysis D: Module Pathway Enrichment\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Load gene sets if provided
if (nchar(args$gene_sets) > 0 && file.exists(args$gene_sets)) {
  if (requireNamespace("yaml", quietly = TRUE)) {
    gs_config <- yaml::yaml.load_file(args$gene_sets)
  } else {
    cat("  WARNING: yaml package not available, skipping pathway enrichment\n")
    gs_config <- NULL
  }
} else {
  cat("  No gene_sets file provided, skipping pathway enrichment\n")
  gs_config <- NULL
}

if (!is.null(gs_config) && nchar(args$id_map) > 0 && file.exists(args$id_map)) {
  # Build symbol → Ensembl map
  id_map <- read.csv(args$id_map, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  eid_col <- grep("ensembl", colnames(id_map), ignore.case = TRUE, value = TRUE)[1]
  sym_col <- grep("symbol|mgi", colnames(id_map), ignore.case = TRUE, value = TRUE)[1]
  sym_to_ens <- setNames(id_map[[eid_col]], tolower(id_map[[sym_col]]))

  # Parse gene sets
  pathway_names <- setdiff(names(gs_config),
                           c("enrichment_databases", "silent_shifter_criteria", "module_detection"))
  background <- colnames(datExpr)  # universe = genes in WGCNA

  enrichment_results <- list()
  for (pw_name in pathway_names) {
    pw <- gs_config[[pw_name]]
    if (!is.list(pw) || is.null(pw$genes)) next
    symbols <- tolower(pw$genes)
    ens_ids <- unique(na.omit(sym_to_ens[symbols]))
    pw_in_bg <- intersect(ens_ids, background)
    if (length(pw_in_bg) < 2) next

    for (mc in unique(moduleColors)) {
      mod_genes <- gene_names[moduleColors == mc]
      overlap <- length(intersect(mod_genes, pw_in_bg))
      # Fisher's exact test (hypergeometric)
      a <- overlap
      b <- length(mod_genes) - overlap
      c_val <- length(pw_in_bg) - overlap
      d <- length(background) - length(mod_genes) - c_val
      ft <- fisher.test(matrix(c(a, b, c_val, max(d, 0)), nrow = 2), alternative = "greater")
      enrichment_results[[length(enrichment_results) + 1]] <- data.frame(
        module = mc,
        pathway = pw_name,
        module_size = length(mod_genes),
        pathway_size_in_bg = length(pw_in_bg),
        overlap = overlap,
        odds_ratio = ft$estimate,
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(enrichment_results) > 0) {
    enrich_df <- do.call(rbind, enrichment_results)
    enrich_df$q_value <- p.adjust(enrich_df$p_value, method = "BH")
    write.csv(enrich_df, file.path(outdir, "module_pathway_enrichment.csv"), row.names = FALSE)

    sig_enrich <- enrich_df[enrich_df$q_value < 0.10 & enrich_df$module != "grey", ]
    if (nrow(sig_enrich) > 0) {
      cat("  Significant enrichments (q < 0.10):\n")
      sig_enrich <- sig_enrich[order(sig_enrich$q_value), ]
      for (i in 1:min(nrow(sig_enrich), 20)) {
        cat(sprintf("    [%s] %s: overlap=%d, OR=%.1f, q=%.4f\n",
                    sig_enrich$module[i], sig_enrich$pathway[i],
                    sig_enrich$overlap[i], sig_enrich$odds_ratio[i],
                    sig_enrich$q_value[i]))
      }
    } else {
      cat("  No pathway enrichments at q < 0.10\n")
    }
  }
} else {
  cat("  Skipping pathway enrichment (missing id_map or gene_sets)\n")
}



summary_info <- list(
  timestamp = Sys.time(),
  n_samples_total = nrow(datExpr),
  n_flt = sum(meta_fg$EnvGroup == "FLT"),
  n_gc = sum(meta_fg$EnvGroup == "GC"),
  n_genes = ncol(datExpr),
  soft_power = chosen_power,
  n_modules = n_modules,
  min_module_size = args$min_module_size,
  merge_cut = args$merge_cut,
  n_pres_perms = args$n_pres_perms,
  network_type = "signed",
  cor_function = "bicor"
)
writeLines(
  paste(names(summary_info), summary_info, sep = ": "),
  file.path(outdir, "wgcna_summary.txt")
)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("WGCNA analysis complete\n")
cat("Outputs in:", outdir, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
