
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# Paths
raw_dir <- "data/processed/deconvolution/raw_ref"
mtx_path <- file.path(raw_dir, "matrix.mtx")
obs_path <- file.path(raw_dir, "barcodes.tsv")
var_path <- file.path(raw_dir, "features.tsv")

message("Loading raw data...")
raw_mat <- readMM(mtx_path)
obs <- data.frame(fread(obs_path, header=TRUE), row.names=1, check.names=FALSE)
var <- data.frame(fread(var_path, header=TRUE), row.names=1, check.names=FALSE)

# Transpose check (Simulating the script's logic)
message("Dimensions before transpose: ", paste(dim(raw_mat), collapse=" x "))
if (ncol(raw_mat) == nrow(var)) {
  message("Transposing...")
  raw_mat <- t(raw_mat)
}
colnames(raw_mat) <- rownames(obs)
rownames(raw_mat) <- rownames(var)

# Redirect output to file
sink("debug_podo_results.txt")

# Target Genes (Ensembl)
# Nphs1, Nphs2, Wt1, Synpo
targets <- c("ENSMUSG00000006649", "ENSMUSG00000026602", "ENSMUSG00000016458", "ENSMUSG00000043079")
names(targets) <- c("Nphs1", "Nphs2", "Wt1", "Synpo")

# Check if targets exist in rownames
found_targets <- intersect(targets, rownames(raw_mat))
cat("\nFound ", length(found_targets), " of ", length(targets), " target genes in matrix rows.\n")
print(found_targets)

# Target Cells
# "Epcam    podocyte" from free_annotation
label_col <- "free_annotation"
if (!label_col %in% colnames(obs)) stop("free_annotation not found in obs")

podo_cells <- rownames(obs)[obs[[label_col]] == "Epcam    podocyte"]
cat("\nFound ", length(podo_cells), " Podocyte cells.\n")

if (length(podo_cells) == 0) stop("No podocytes found!")

# Check expression of targets in podocyte cells
cat("\n--- Expression in Podocytes (Raw Counts) ---\n")
podo_mat <- raw_mat[found_targets, podo_cells, drop=FALSE]

for (g_id in found_targets) {
  g_name <- names(targets)[targets == g_id]
  vals <- podo_mat[g_id, ]
  
  cat(sprintf("Gene %s (%s):\n", g_name, g_id))
  cat(sprintf("  Non-zero cells: %d / %d (%.1f%%)\n", sum(vals > 0), length(vals), 100 * sum(vals > 0)/length(vals)))
  cat(sprintf("  Mean count: %.4f\n", mean(vals)))
  cat(sprintf("  Max count:  %d\n", max(vals)))
}

# Check random other cells to ensure it's not simply swapped
cat("\n--- Expression in random NON-podocytes (Control) ---\n")
non_podo_cells <- rownames(obs)[obs[[label_col]] != "Epcam    podocyte"]
if (length(non_podo_cells) > 1000) non_podo_cells <- sample(non_podo_cells, 1000)

ctrl_mat <- raw_mat[found_targets, non_podo_cells, drop=FALSE]

for (g_id in found_targets) {
  g_name <- names(targets)[targets == g_id]
  vals <- ctrl_mat[g_id, ]
  
  cat(sprintf("Gene %s (%s):\n", g_name, g_id))
  cat(sprintf("  Non-zero cells: %d / %d (%.1f%%)\n", sum(vals > 0), length(vals), 100 * sum(vals > 0)/length(vals)))
  cat(sprintf("  Mean count: %.4f\n", mean(vals)))
}
sink()
