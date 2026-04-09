# Install all required R packages for RRRM-2 Pipeline
# Run this script once: Rscript install_r_packages.R

cat("=== Installing R packages for RRRM-2 Pipeline ===\n\n")

# 1) Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 2) Bioconductor packages
bioc_packages <- c(
    "SingleCellExperiment",
    "Biobase",
    "edgeR",
    "DESeq2",
    "AnnotationDbi",
    "org.Mm.eg.db",
    "zellkonverter",
    "sva"
)

cat("Installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("  Installing: ", pkg, "\n"))
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
        cat(paste0("  Already installed: ", pkg, "\n"))
    }
}

# 3) CRAN packages
cran_packages <- c(
    "data.table",
    "Matrix",
    "ggplot2"
)

cat("\nInstalling CRAN packages...\n")
for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("  Installing: ", pkg, "\n"))
        install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
        cat(paste0("  Already installed: ", pkg, "\n"))
    }
}

# 4) MuSiC package (GitHub)
if (!require("MuSiC", quietly = TRUE)) {
    cat("\nInstalling MuSiC from GitHub...\n")
    if (!require("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github("xuranw/MuSiC")
} else {
    cat("\nAlready installed: MuSiC\n")
}

# 5) Optional: ashr for LFC shrinkage in differential expression
if (!require("ashr", quietly = TRUE)) {
    cat("\nInstalling optional ashr package...\n")
    install.packages("ashr", repos = "https://cloud.r-project.org")
} else {
    cat("\nAlready installed: ashr (optional)\n")
}

cat("\n=== Verifying installations ===\n")
all_packages <- c(bioc_packages, cran_packages, "MuSiC")
success <- sapply(all_packages, function(pkg) {
    loaded <- require(pkg, character.only = TRUE, quietly = TRUE)
    cat(paste0("  ", pkg, ": ", ifelse(loaded, "OK", "FAILED"), "\n"))
    loaded
})

if (all(success)) {
    cat("\n✓ All required packages installed successfully!\n")
} else {
    cat("\n✗ Some packages failed to install. See above for details.\n")
    cat("Failed packages:", paste(names(success)[!success], collapse = ", "), "\n")
}
