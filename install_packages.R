# COLABS-IRTLab Project - Package Installation Script
# This script installs all required packages for the project
# Run this script to set up your R environment

# --- Core CRAN Packages ---
cran_packages <- c(
  # Data manipulation and analysis
  "dplyr",
  "tidyr", 
  "tidyverse",
  "stringr",
  
  # Visualization
  "ggplot2",
  "gplots",
  "ggpattern",
  "RColorBrewer",
  "colorspace",
  "viridis",
  "ComplexHeatmap",
  "circlize",
  "dendextend",
  "dichromat",
  "Polychrome",
  
  # Date/time handling
  "lubridate",
  
  # Data import/export
  "readxl",
  
  # Statistical modeling
  "growthcurver",
  
  # Parallel processing
  "doSNOW",
  "foreach",
  "progress",
  
  # Utilities
  "gtools",
  "slider", 
  "grid",
  
  # Additional packages from Quarto files
  "forcats",
  "glue",
  "scales"
)

# --- Bioconductor Packages ---
bioc_packages <- c(
  # Sequence analysis
  "Biostrings",
  "ape", 
  "seqinr",
  "pwalign",
  
  # DADA2 pipeline and sequence analysis (from Quarto files)
  "dada2",
  "DECIPHER"
)

# Function to install packages if not already installed
install_if_missing <- function(packages, repo = "CRAN") {
  if (repo == "CRAN") {
    missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
    if (length(missing_packages) > 0) {
      cat("Installing CRAN packages:", paste(missing_packages, collapse = ", "), "\n")
      install.packages(missing_packages, dependencies = TRUE)
    } else {
      cat("All CRAN packages already installed.\n")
    }
  } else if (repo == "Bioconductor") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
    if (length(missing_packages) > 0) {
      cat("Installing Bioconductor packages:", paste(missing_packages, collapse = ", "), "\n")
      BiocManager::install(missing_packages, dependencies = TRUE)
    } else {
      cat("All Bioconductor packages already installed.\n")
    }
  }
}

# Install packages
cat("=== Installing COLABS-IRTLab Project Dependencies ===\n\n")

cat("1. Installing CRAN packages...\n")
install_if_missing(cran_packages, "CRAN")

cat("\n2. Installing Bioconductor packages...\n")
install_if_missing(bioc_packages, "Bioconductor")

cat("\n=== Installation Complete! ===\n")
cat("You can now run the R scripts in the project.\n\n")

# Verify installation by loading key packages
cat("Verifying installation of key packages...\n")
key_packages <- c("dplyr", "ggplot2", "Biostrings", "ape")
for (pkg in key_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗", pkg, "failed to load\n")
  }
}
