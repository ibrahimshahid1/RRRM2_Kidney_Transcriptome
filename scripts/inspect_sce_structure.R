
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
})

h5ad_path <- "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"

message("Loading H5AD: ", h5ad_path)
sce <- readH5AD(h5ad_path)

message("\n--- Assays ---")
print(assayNames(sce))

message("\n--- Alt Experiments ---")
print(altExpNames(sce))

message("\n--- Reduced Dims ---")
print(reducedDimNames(sce))

# Check samples of X
message("\n--- X (first 5x5) ---")
print(assay(sce, "X")[1:5, 1:5])

# Check if 'raw' or 'counts' exists and inspect
if ("raw" %in% assayNames(sce)) {
    message("\n--- Raw Assay (first 5x5) ---")
    print(assay(sce, "raw")[1:5, 1:5])
}

if ("counts" %in% assayNames(sce)) {
    message("\n--- Counts Assay (first 5x5) ---")
    print(assay(sce, "counts")[1:5, 1:5])
}

# Check raw in altExp if presumed there
# Often scanpy raw accounts are distinct
