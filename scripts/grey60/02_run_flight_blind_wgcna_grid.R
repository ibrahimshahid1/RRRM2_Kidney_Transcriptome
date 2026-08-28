#!/usr/bin/env Rscript
# Rebuild Grey60 candidates using only BSL+VIV reference animals.

suppressPackageStartupMessages({
  library(WGCNA)
  library(optparse)
})

disableWGCNAThreads()
options(stringsAsFactors = FALSE)

option_list <- list(
  make_option("--rtech", type = "character", required = TRUE),
  make_option("--meta", type = "character", required = TRUE),
  make_option("--outdir", type = "character", required = TRUE),
  make_option("--power", type = "integer", default = 7),
  make_option("--seed", type = "integer", default = 7716029)
)
args <- parse_args(OptionParser(option_list = option_list))
set.seed(args$seed)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
grid_dir <- file.path(args$outdir, "grid_assignments")
dir.create(grid_dir, recursive = TRUE, showWarnings = FALSE)

rtech <- read.csv(
  gzfile(args$rtech), sep = "\t", row.names = 1, check.names = FALSE
)
meta <- read.csv(gzfile(args$meta), sep = "\t", check.names = FALSE)
sample_col <- NA
for (candidate in c("Sample Name (raw_counts_colname)", "Sample Name", "sample")) {
  if (candidate %in% colnames(meta)) {
    sample_col <- candidate
    break
  }
}
if (is.na(sample_col)) stop("Could not resolve metadata sample column")
rownames(meta) <- meta[[sample_col]]
meta$EnvGroup[meta$EnvGroup == "HGC"] <- "GC"
meta$EnvGroup[meta$EnvGroup == "VGC"] <- "VIV"

common <- intersect(colnames(rtech), rownames(meta))
rtech <- rtech[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
reference <- meta$EnvGroup %in% c("BSL", "VIV")
if (sum(reference) != 40) {
  stop(sprintf("Expected 40 BSL+VIV reference animals; found %d", sum(reference)))
}
ref_expr <- rtech[, reference, drop = FALSE]

max_genes_values <- c(2500, 5000, 10000)
min_module_values <- c(20, 30, 50)
merge_values <- c(0.15, 0.25, 0.35)
summary_rows <- list()
row_i <- 1

for (max_genes in max_genes_values) {
  gene_var <- sort(apply(ref_expr, 1, var, na.rm = TRUE), decreasing = TRUE)
  selected <- names(gene_var)[seq_len(min(max_genes, length(gene_var)))]
  datExpr <- as.data.frame(t(ref_expr[selected, , drop = FALSE]))
  qc <- goodSamplesGenes(datExpr, verbose = 0)
  if (!qc$allOK) datExpr <- datExpr[qc$goodSamples, qc$goodGenes, drop = FALSE]

  for (min_module in min_module_values) {
    for (merge_cut in merge_values) {
      variant <- sprintf(
        "g%05d_m%02d_c%02d",
        max_genes, min_module, as.integer(round(merge_cut * 100))
      )
      output <- file.path(grid_dir, paste0(variant, ".tsv.gz"))
      started <- Sys.time()
      status <- "completed"
      error_message <- ""
      module_count <- NA_integer_

      if (!file.exists(output)) {
        result <- tryCatch({
          net <- blockwiseModules(
            datExpr,
            power = args$power,
            networkType = "signed",
            TOMType = "signed",
            corType = "bicor",
            maxPOutliers = 0.05,
            minModuleSize = min_module,
            mergeCutHeight = merge_cut,
            numericLabels = TRUE,
            pamRespectsDendro = FALSE,
            saveTOMs = FALSE,
            maxBlockSize = 5000,
            randomSeed = args$seed,
            verbose = 0
          )
          colors <- labels2colors(net$colors)
          module_count <- length(setdiff(unique(colors), "grey"))
          assignments <- data.frame(
            gene = colnames(datExpr),
            module_num = net$colors,
            module_color = colors,
            stringsAsFactors = FALSE
          )
          con <- gzfile(output, "wt")
          write.table(
            assignments, con, sep = "\t", row.names = FALSE, quote = FALSE
          )
          close(con)
          TRUE
        }, error = function(e) {
          status <<- "failed"
          error_message <<- conditionMessage(e)
          FALSE
        })
      } else {
        status <- "reused"
        existing <- read.csv(gzfile(output), sep = "\t")
        module_count <- length(setdiff(unique(existing$module_color), "grey"))
      }

      elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      summary_rows[[row_i]] <- data.frame(
        variant = variant,
        max_genes = max_genes,
        min_module_size = min_module,
        merge_cut = merge_cut,
        power = args$power,
        n_reference_samples = nrow(datExpr),
        n_genes = ncol(datExpr),
        n_modules = module_count,
        status = status,
        elapsed_seconds = elapsed,
        output = output,
        error = error_message,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1
      write.table(
        do.call(rbind, summary_rows),
        file.path(args$outdir, "grid_run_summary.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
      )
      cat(sprintf(
        "[%s] %s: modules=%s elapsed=%.1fs\n",
        format(Sys.time(), "%H:%M:%S"), variant, module_count, elapsed
      ))
    }
  }
}

cat("Flight-blind WGCNA grid complete:", args$outdir, "\n")
