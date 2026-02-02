# Install limma from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("limma", ask = FALSE, update = FALSE)
cat("limma installation complete!\n")
