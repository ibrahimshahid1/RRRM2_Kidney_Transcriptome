# Marker-based cluster annotation for Chen atlas
suppressPackageStartupMessages(library(Matrix))

message("Loading Chen atlas...")
load("data/external/single_cell_atlases/chen_atlas/GSE150338_Seurat_IntegrateData.RData")

obj <- merge.data.integrated
counts <- obj@assays$RNA@counts
meta <- obj@meta.data
genes <- rownames(counts)

# Markers to check
markers <- c(
    "Slc12a1", "Umod", "Slc12a3", "Pvalb", "Trpm6",
    "Calb1", "Aqp2", "Slc26a4", "Slc34a1",
    "Egf", "Pth1r", "Clcnkb", "Scnn1g", "Wnt4"
)

found <- markers[markers %in% genes]
missing <- markers[!markers %in% genes]
message("Found markers: ", paste(found, collapse = ", "))
message("Missing: ", paste(missing, collapse = ", "))

# Compute mean expression per cluster for each marker
clusters <- meta$seurat_clusters
cluster_levels <- sort(unique(clusters))

result <- matrix(0, nrow = length(found), ncol = length(cluster_levels))
rownames(result) <- found
colnames(result) <- cluster_levels

for (g in found) {
    for (cl in cluster_levels) {
        cells <- which(clusters == cl)
        result[g, as.character(cl)] <- mean(counts[g, cells])
    }
}

message("\n--- Mean expression of markers per cluster ---")
print(round(result, 2))

# Also show detection rate (% of cells expressing)
detect <- matrix(0, nrow = length(found), ncol = length(cluster_levels))
rownames(detect) <- found
colnames(detect) <- cluster_levels

for (g in found) {
    for (cl in cluster_levels) {
        cells <- which(clusters == cl)
        detect[g, as.character(cl)] <- round(100 * mean(counts[g, cells] > 0), 1)
    }
}

message("\n--- Detection rate (%) of markers per cluster ---")
print(detect)

# Print cluster sizes for reference
message("\n--- Cluster sizes ---")
print(table(clusters))

message("\nDone.")
