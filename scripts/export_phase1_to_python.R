# scripts/export_phase1_to_python.R
# Export Phase 1 residualized expression data to Python-friendly formats

obj <- readRDS("data/processed/phase1_residuals/phase1_Rtech.rds")

Rtech <- obj$Rtech # genes x samples
meta <- obj$meta

outdir <- "data/processed/phase1_residuals"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Expression
gz1 <- gzfile(file.path(outdir, "Rtech.tsv.gz"), "wt")
write.table(Rtech, gz1, sep = "\t", quote = FALSE, col.names = NA)
close(gz1)

# Metadata (already aligned to sample order)
gz2 <- gzfile(file.path(outdir, "meta_phase1.tsv.gz"), "wt")
write.table(meta, gz2, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz2)

cat(
    "Wrote:\n",
    file.path(outdir, "Rtech.tsv.gz"), "\n",
    file.path(outdir, "meta_phase1.tsv.gz"), "\n"
)
