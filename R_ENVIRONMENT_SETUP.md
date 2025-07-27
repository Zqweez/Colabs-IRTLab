# COLABS-IRTLab Project Environment Setup
# Created: July 27, 2025

This document describes the R package dependencies for the COLABS-IRTLab project and provides multiple ways to reproduce the computational environment.

## Quick Start

### Option 1: Automatic Installation (Recommended)
```r
source("install_packages.R")
```

### Option 2: Reproducible Environment with renv (Best for Collaboration)
```r
source("setup_renv.R")
```

### Option 3: Manual Installation
```r
# Install required packages manually (see requirements.txt for full list)
install.packages(c("tidyverse", "ggplot2", "ape", "seqinr", ...))
BiocManager::install(c("Biostrings", "pwalign"))
```

## Package Dependencies by Analysis Type

### NGS Analysis (`scripts/ngs/`)
- **Core packages**: `ape`, `seqinr`, `Biostrings`, `pwalign`
- **Data manipulation**: `dplyr`, `tidyr`, `stringr`
- **Visualization**: `ggplot2`, `dendextend`
- **Parallel processing**: `doSNOW`, `foreach`, `progress`

### Growth Curve Analysis (`scripts/growth/`)
- **Core packages**: `tidyverse`, `growthcurver`
- **Visualization**: `ggplot2`, `gplots`
- **Data handling**: `lubridate`, `slider`

### Sanger Sequencing Analysis (`scripts/sanger/`)
- **Sequence analysis**: `ape`, `seqinr`, `Biostrings`, `pwalign`
- **Clustering**: `dendextend`, `ComplexHeatmap`
- **Visualization**: `ggplot2`, `viridis`, `RColorBrewer`, `circlize`
- **Parallel processing**: `doSNOW`, `foreach`

### GCMS Analysis (`scripts/gcms/`)
- **Data manipulation**: `dplyr`, `readxl`, `gtools`
- **Visualization**: `ggplot2`, `ggpattern`
- **Date handling**: `lubridate`

## Environment Management Options

### 1. renv (Recommended for Reproducibility)
`renv` creates project-specific package libraries and lock files:

**Setup:**
```r
source("setup_renv.R")
```

**Daily workflow:**
```r
renv::snapshot()  # Save current state
renv::restore()   # Restore from lock file
renv::status()    # Check for changes
```

**Benefits:**
- Exact package versions recorded
- Isolated from system R library
- Perfect for collaboration
- Version control friendly

### 2. Traditional Package Installation
Install packages globally in your R installation:

```r
source("install_packages.R")
```

**Benefits:**
- Simple and familiar
- Packages available across all projects
- No setup required

## System Requirements

- **R version**: >= 4.0.0 (recommended: >= 4.2.0)
- **Bioconductor**: Will be automatically installed if needed
- **RAM**: >= 8GB recommended for sequence alignment tasks
- **CPU cores**: Multiple cores recommended for parallel processing (scripts use `parallel::detectCores() - 2`)

## Troubleshooting

### Common Issues:

1. **Bioconductor packages not installing**:
   ```r
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("Biostrings")
   ```

2. **Parallel processing errors**:
   - Reduce number of cores in scripts if you encounter memory issues
   - Modify `num_cores <- max(1L, parallel::detectCores() - 2L)` to use fewer cores

3. **Package compilation issues on macOS**:
   - Install Xcode command line tools: `xcode-select --install`
   - Install gfortran if needed

4. **Memory issues with large datasets**:
   - Close other applications
   - Increase virtual memory/swap space
   - Consider running analyses on subsets of data

## File Structure

```
├── install_packages.R      # Automated package installation
├── setup_renv.R           # renv environment setup
├── requirements.txt       # Human-readable package list
├── renv.lock             # Auto-generated lock file (if using renv)
└── renv/                 # renv directory (if using renv)
```

## Sharing Your Environment

### With renv:
1. Commit `renv.lock` to version control
2. Collaborators run: `renv::restore()`

### Without renv:
1. Share `requirements.txt` and `install_packages.R`
2. Collaborators run: `source("install_packages.R")`

## Package Sources

- **CRAN packages**: Standard R package repository
- **Bioconductor packages**: Specialized bioinformatics packages
- All packages are from official repositories (no GitHub/custom sources)

For questions about specific package functionality, refer to:
- Individual package documentation: `?packagename`
- Bioconductor documentation: https://bioconductor.org/
- CRAN documentation: https://cran.r-project.org/
