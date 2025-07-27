# COLABS-IRTLab Project - renv Setup Script
# This script initializes a reproducible R environment using renv
# 
# renv is the modern R package management system that creates project-specific
# libraries and lock files for reproducibility across different machines/users
#
# Usage:
# 1. Run this script once: source("setup_renv.R")
# 2. Install packages normally, renv will track them
# 3. Use renv::snapshot() to save current state
# 4. Share renv.lock file with collaborators
# 5. Collaborators run renv::restore() to get exact same environment

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv for this project
if (!file.exists("renv.lock")) {
  cat("Initializing renv for COLABS-IRTLab project...\n")
  renv::init()
} else {
  cat("renv already initialized. Use renv::restore() to install dependencies.\n")
}

# Install all required packages
cat("Installing project dependencies...\n")

# CRAN packages
cran_packages <- c(
  "dplyr", "tidyr", "tidyverse", "stringr",
  "ggplot2", "gplots", "ggpattern", "RColorBrewer", 
  "colorspace", "viridis", "ComplexHeatmap", "circlize",
  "dendextend", "dichromat", "Polychrome", "lubridate",
  "readxl", "growthcurver", "doSNOW", "foreach", 
  "progress", "gtools", "slider", "grid", "ape", "seqinr",
  "forcats", "glue", "scales"
)

# Install CRAN packages
install.packages(cran_packages)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("Biostrings", "pwalign", "dada2", "DECIPHER"))

# Take a snapshot of the current state
cat("Creating renv snapshot...\n")
renv::snapshot()

cat("\n=== renv Setup Complete! ===\n")
cat("Your R environment is now reproducible.\n")
cat("Key commands:\n")
cat("- renv::snapshot(): Save current package state\n")
cat("- renv::restore(): Restore packages from lock file\n")
cat("- renv::status(): Check for changes\n")
cat("- renv::update(): Update packages\n\n")

cat("To share this environment:\n")
cat("1. Share the renv.lock file with collaborators\n")
cat("2. They should run: renv::restore()\n")
