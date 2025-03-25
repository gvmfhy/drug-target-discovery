#!/usr/bin/env Rscript

# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to install a package if it's not already installed
install_if_missing <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing package:", package_name))
    tryCatch({
      if (package_name %in% rownames(available.packages())) {
        install.packages(package_name)
      } else {
        BiocManager::install(package_name, update = FALSE, ask = FALSE)
      }
    }, error = function(e) {
      message(paste("Error installing", package_name, ":", e$message))
    })
  } else {
    message(paste("Package", package_name, "is already installed"))
  }
}

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# First install base R dependencies
base_packages <- c(
  "jsonlite",
  "curl",
  "httr",
  "xml2",
  "Matrix",
  "statmod"
)

message("Installing base R dependencies...")
for (package in base_packages) {
  install_if_missing(package)
}

# List of required Bioconductor packages
bioc_packages <- c(
  "BiocGenerics",
  "Biobase",
  "S4Vectors",
  "IRanges",
  "GenomeInfoDb",
  "AnnotationDbi",
  "Biostrings",
  "GenomicRanges",
  "XVector",
  "zlibbioc",
  "hgu133plus2.db",  # For Affymetrix HG-U133 Plus 2.0 array annotation
  "limma"            # For differential expression analysis
)

message("\nInstalling Bioconductor packages...")
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# Verify installations
message("\nVerifying installations...")
for (package in c(base_packages, bioc_packages)) {
  if (require(package, character.only = TRUE, quietly = TRUE)) {
    message(paste(package, "successfully installed"))
  } else {
    message(paste("WARNING:", package, "installation may have failed"))
  }
} 