#!/usr/bin/env Rscript
# GC-only reference WGCNA → FLT preservation test
# Item 5: Verify preservation with proper GC-reference modules

suppressPackageStartupMessages(library(WGCNA))
disableWGCNAThreads()
options(stringsAsFactors = FALSE)

cat("============================================================\n")
cat("GC-Only Reference WGCNA + Preservation Test\n")
cat("============================================================\n\n")

args <- commandArgs(trailingOnly = TRUE)
rtech_path <- args[1]
meta_path  <- args[2]
outdir     <- args[3]
n_perms    <- ifelse(length(args) >= 4, as.integer(args[4]), 200)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# Load data
rtech <- read.csv(gzfile(rtech_path), sep = "\t", row.names = 1, check.names = FALSE)
meta <- read.csv(gzfile(meta_path), sep = "\t", check.names = FALSE)
sample_col <- NA
for (col in c("Sample Name (raw_counts_colname)", "Sample Name", "sample")) {
  if (col %in% colnames(meta)) { sample_col <- col; break }
}
if (is.na(sample_col)) sample_col <- colnames(meta)[1]
rownames(meta) <- meta[[sample_col]]
meta$EnvGroup[meta$EnvGroup == "HGC"] <- "GC"
meta$Age[meta$Age %in% c("Young", "young", "Yng")] <- "YNG"
meta$Age[meta$Age %in% c("old")] <- "OLD"
meta$Arm[meta$Arm %in% c("ISS", "ISST", "ISS_T")] <- "ISS-T"

common <- intersect(colnames(rtech), rownames(meta))
rtech <- rtech[, common]
meta <- meta[common, ]

# FLT + GC only
keep <- meta$EnvGroup %in% c("FLT", "GC")
rtech_fg <- rtech[, keep]
meta_fg <- meta[keep, ]

# Top 5000 by variance
gene_var <- sort(apply(rtech_fg, 1, var, na.rm = TRUE), decreasing = TRUE)
top_genes <- names(gene_var)[1:min(5000, length(gene_var))]

gc_mask  <- meta_fg$EnvGroup == "GC"
flt_mask <- meta_fg$EnvGroup == "FLT"

datExpr_GC  <- as.data.frame(t(rtech_fg[top_genes, gc_mask]))
datExpr_FLT <- as.data.frame(t(rtech_fg[top_genes, flt_mask]))

cat("GC samples:", nrow(datExpr_GC), "\n")
cat("FLT samples:", nrow(datExpr_FLT), "\n")
cat("Genes:", ncol(datExpr_GC), "\n\n")

# Build WGCNA on GC-only
cat("Building GC-only network...\n")
sft <- pickSoftThreshold(datExpr_GC, powerVector = c(1:20, seq(22, 30, by = 2)),
                         networkType = "signed", corFnc = "bicor",
                         corOptions = list(maxPOutliers = 0.05), verbose = 0)

fit_vals <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
above <- which(fit_vals > 0.80)
power <- if (length(above) > 0) sft$fitIndices[above[1], 1] else 12
cat("  GC-only soft threshold power:", power, "\n")

net_gc <- blockwiseModules(datExpr_GC, power = power, networkType = "signed",
                           TOMType = "signed", corType = "bicor",
                           maxPOutliers = 0.05, minModuleSize = 30,
                           mergeCutHeight = 0.25, numericLabels = TRUE,
                           pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0)

gc_colors <- labels2colors(net_gc$colors)
cat("  GC-only modules:", length(unique(gc_colors)) - ("grey" %in% gc_colors), "\n")
cat("  Module sizes:\n")
print(sort(table(gc_colors), decreasing = TRUE))

# Save GC-only module assignments
gc_mod_df <- data.frame(gene = colnames(datExpr_GC),
                        module_num_gc = net_gc$colors,
                        module_color_gc = gc_colors)
write.csv(gc_mod_df, file.path(outdir, "gc_only_module_assignments.csv"), row.names = FALSE)

# Preservation: GC reference → FLT test
cat("\nRunning GC→FLT preservation (", n_perms, "permutations)...\n")

multiExpr <- list(
  GC  = list(data = datExpr_GC),
  FLT = list(data = datExpr_FLT)
)
multiColor <- list(GC = gc_colors)

pres <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                           nPermutations = n_perms, networkType = "signed",
                           verbose = 0)

pres_stats <- pres$preservation$Z$ref.GC$inColumnsAlsoPresentIn.FLT
pres_obs <- pres$preservation$observed$ref.GC$inColumnsAlsoPresentIn.FLT

# Safe column access
get_col <- function(df, patterns) {
  for (p in patterns) {
    m <- grep(p, colnames(df), value = TRUE, ignore.case = TRUE)
    if (length(m) > 0) return(df[, m[1]])
  }
  return(rep(NA, nrow(df)))
}

pres_df <- data.frame(
  module_color_gc = rownames(pres_stats),
  moduleSize = get_col(pres_obs, c("moduleSize")),
  Zsummary = get_col(pres_stats, c("Zsummary.pres", "Zsummary")),
  Zdensity = get_col(pres_stats, c("Zdensity.pres", "Zdensity")),
  Zconnectivity = get_col(pres_stats, c("Zconnectivity.pres", "Zconnectivity"))
)
pres_df <- pres_df[!pres_df$module_color_gc %in% c("gold"), ]
pres_df$preservation <- ifelse(is.na(pres_df$Zsummary), "unknown",
                        ifelse(pres_df$Zsummary > 10, "strong",
                        ifelse(pres_df$Zsummary > 2, "moderate", "none")))

write.csv(pres_df, file.path(outdir, "gc_reference_preservation.csv"), row.names = FALSE)

cat("\nGC-reference preservation results:\n")
pres_df <- pres_df[order(-pres_df$moduleSize), ]
for (i in 1:nrow(pres_df)) {
  cat(sprintf("  %-15s (n=%3d): Zsummary=%6.2f  [%s]\n",
              pres_df$module_color_gc[i],
              ifelse(is.na(pres_df$moduleSize[i]), 0, pres_df$moduleSize[i]),
              ifelse(is.na(pres_df$Zsummary[i]), 0, pres_df$Zsummary[i]),
              pres_df$preservation[i]))
}

# Cross-reference: which combined-WGCNA module genes fall into which GC-only modules?
# Load combined module assignments
combined_path <- file.path(dirname(outdir), "module_assignments.csv")
if (file.exists(combined_path)) {
  combined <- read.csv(combined_path)
  merged <- merge(combined[, c("gene", "module_color")],
                  gc_mod_df[, c("gene", "module_color_gc")], by = "gene")
  
  cat("\nCross-reference: combined → GC-only module overlap:\n")
  for (mc in c("grey60", "red", "midnightblue", "blue")) {
    sub <- merged[merged$module_color == mc, ]
    if (nrow(sub) > 0) {
      gc_dist <- sort(table(sub$module_color_gc), decreasing = TRUE)
      top3 <- paste(paste0(names(gc_dist)[1:min(3, length(gc_dist))], "(",
                           gc_dist[1:min(3, length(gc_dist))], ")"), collapse = ", ")
      cat(sprintf("  %-15s → GC-only: %s\n", mc, top3))
    }
  }
}

cat("\nOutputs in:", outdir, "\n")
