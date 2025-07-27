# COLABS-IRTLab: Microbial Analysis Pipelines

This repository is associated with the work done in the IRTLab at the Laboratory of Microbial Genetics and Evolution, Graduate School of Life Sciences, Tohoku University.

This repository contains comprehensive scripts, data structures, and workflows for analyzing microbial data generated from GC-MS, Sanger sequencing, next-generation sequencing (NGS), and growth curve experiments. The project characterizes microbial isolates, quantifies phenanthrene degradation, performs taxonomic identification, and clusters microbial communities based on sequence similarity and phylogenetic relationships.

------------------------------------------------------------------------

## Project Structure

```
project-root/
├── data/                    # Raw and processed input data
│   ├── gcms_raw/           # GC-MS chromatogram data (.txt, .xlsx, .pdf)
│   ├── growth/             # Growth curve data (.dt files, .csv, .txt)
│   ├── phe_growth/         # Phenanthrene growth experiments
│   ├── ngs/                # NGS sequencing data
│   │   ├── raw_reads/      # Raw FASTQ files
│   │   ├── filtered/       # Filtered reads
│   │   └── fastqc/         # Quality control reports
│   └── sanger/             # Sanger sequencing data and results
│       ├── 00_ab1_raw_reads/           # Raw .ab1 files
│       ├── 01_ab1_trimmed_reads/       # Trimmed sequences
│       ├── 02_ab1_aligned_read_pairs/  # Aligned read pairs
│       ├── 03_ab1_assembled_sequences/ # Consensus sequences
│       └── 04_taxa_files/              # Taxonomic assignments
├── databases/              # Reference databases
│   ├── ASVs/              # ASV reference database
│   ├── ezbiocloud/        # EZBioCloud 16S database
│   ├── Isolates/          # Isolate consensus sequences
│   └── SILVA/             # SILVA 16S database
├── results/                # Analysis outputs
│   ├── gcms_plots/        # GC-MS visualization outputs
│   ├── growth_models/     # Growth curve analysis results
│   ├── growth_phe/        # Phenanthrene growth results
│   ├── ngs/               # NGS analysis outputs
│   └── sanger_plots/      # Sanger sequencing visualizations
├── scripts/                # Analysis scripts
│   ├── gcms/              # GC-MS analysis (R)
│   ├── growth/            # Growth curve analysis (R)
│   ├── ngs/               # NGS pipelines (Python, R, Quarto)
│   └── sanger/            # Sanger sequencing analysis (R)
├── install_packages.R      # Automated R package installation
├── setup_renv.R           # Reproducible environment setup
├── R_ENVIRONMENT_SETUP.md  # Detailed environment guide
├── README.md              # This file
└── 00 Colabs-repo.Rproj   # RStudio project file
```

------------------------------------------------------------------------

## Analysis Workflows

### 1. GC-MS Analysis (Phenanthrene Degradation)

**Goal**: Quantify phenanthrene degradation using GC-MS chromatography.

**Data**:
- `data/gcms_raw/L-tubes.xlsx`: Processed phenanthrene/decane ratios
- `data/gcms_raw/naphtol-calibration.xlsx`: 1-Naphthol calibration data
- `data/gcms_raw/*.txt`: Raw chromatogram outputs
- `data/gcms_raw/*.pdf`: Chromatogram visualizations

**Scripts**:
- `scripts/gcms/GCMS-phe-nap.R`: Analyzes phenanthrene degradation across samples and dates
- `scripts/gcms/GCMS-nap-calibration.R`: Naphthol calibration curve analysis

**Output**: `results/gcms_plots/` - Bar plots with error bars showing degradation over time

---

### 2. Growth Curve Analysis

**Goal**: Analyze microbial growth using OD600 measurements from biophotorecorder.

**Data**:
- `data/growth/*.dt*`: Raw biophotorecorder data files
- `data/growth/*.csv`: Converted growth data with headers
- `data/phe_growth/*.dt*`: Growth with phenanthrene supplementation

**Scripts**:
- `scripts/growth/plot-growth.R`: Basic growth curve visualization and statistical modeling
- `scripts/growth/plot-growth-phe.R`: Phenanthrene-supplemented growth analysis

**Features**:
- Automatic time correction for midnight rollovers
- Sliding window smoothing
- Growth parameter estimation (μ_max, lag time, generation time)

**Output**: `results/growth_models/`, `results/growth_phe/` - Growth curves and statistical summaries

---

### 3. Sanger Sequencing Analysis

**Goal**: Process .ab1 files, generate consensus sequences, assign taxonomy, and perform phylogenetic clustering.

**Data**:
- `data/sanger/00_ab1_raw_reads/`: Raw .ab1 chromatogram files
- `data/sanger/formatted_consensus.fasta`: Final consensus sequences
- `data/sanger/blast_top_hit_*.tsv`: Taxonomic assignments
- `data/sanger/taxa-genus-species-pid.csv`: Taxa with percent identity

**Scripts**:
- `scripts/sanger/16s_cluster_kmers.R`: K-mer based sequence clustering and PCA
- `scripts/sanger/16s_cluster_alignment.R`: Pairwise alignment clustergrams

**Features**:
- Parallel pairwise sequence alignment
- Weighted percent identity calculations
- Hierarchical clustering with multiple linkage methods
- Complex heatmap visualization with taxonomic annotations

**Output**: `results/sanger_plots/` - Dendrograms, heatmaps, PCA plots

---

### 4. NGS Analysis (16S rRNA Amplicon Sequencing)

**Goal**: Process Illumina paired-end reads using DADA2 pipeline for ASV identification and community analysis.

**Data**:
- `data/ngs/raw_reads/`: Raw FASTQ files organized by sample type
- `data/ngs/filtered/`: Quality-filtered reads
- `data/ngs/fastqc/`: Quality control reports

**Scripts**:
- `scripts/ngs/ASV.qmd`: Complete DADA2 pipeline (QC → filtering → denoising → ASV table)
- `scripts/ngs/analyze_asvs.qmd`: ASV community analysis, clustering, and visualization
- `scripts/ngs/run_blast_ngs.py`: Taxonomic assignment using BLAST against multiple databases
- `scripts/ngs/16S-comparison-KS-vs-original.R`: Compare ASVs with Sanger consensus sequences

**Features**:
- Automated quality assessment and filtering
- Error learning and denoising
- Chimera removal
- ASV clustering based on sequence similarity
- Taxonomic assignment with multiple databases (SILVA, EZBioCloud)
- Community composition analysis with diversity metrics

**Output**: `results/ngs/` - ASV tables, taxonomic assignments, community plots

---

### 5. Comparative Analysis

**Goal**: Integrate results across different sequencing methods and experimental conditions.

**Scripts**:
- `scripts/ngs/16S-comparison-KS-vs-original.R`: Detailed comparison between Sanger and NGS results
- Cross-reference isolate identities with growth and degradation capabilities

---

## Path Handling

All scripts assume you are running from the **repository root** and use **relative paths**. 

When working in RStudio, ensure your working directory is set to the project root:
```r
getwd()  # Should show the repository root
# If not, use: setwd("path/to/COLABS-IRTLab")
```

---

## Key Features

- **Comprehensive Pipeline**: Complete workflows from raw data to figures
- **Parallel Processing**: Sequence alignment using multiple CPU cores
- **Reproducible Environment**: `renv` integration for consistent package versions
- **Quality Control**: Automated QC checks and filtering of ab1 and NGS data
- **Multiple Databases**: Support for SILVA, EZBioCloud, and custom databases
- **Flexible Analysis**: Modular scripts that can be run independently or as part of larger workflows

---

## Contact
Laboratory of Microbial Genetics and Evolution  
Graduate School of Life Sciences  
Tohoku University
Japan

Life Sciences
Chalmers University of Technology
Sweden

For questions about the code, please open an issue on GitHub. Or contact me at 
kasperst@chalmers.se

------------------------------------------------------------------------

## Getting Started

### Quick Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/Zqweez/COLABS-IRTLab.git
   cd COLABS-IRTLab
   ```

2. **Set up R environment** (choose one option):
   
   **Option A: Automatic installation**
   ```r
   source("install_packages.R")
   ```
   
   **Option B: Reproducible environment with renv**
   ```r
   source("setup_renv.R")
   ```

3. **For detailed environment setup**, see `R_ENVIRONMENT_SETUP.md`

### RStudio Memory Management (macOS Users)

If you're using RStudio on macOS and encounter memory issues during sequence analysis, restrict RAM usage by running this in Terminal before opening RStudio:

```bash
export R_MAX_VSIZE=16Gi
open -a RStudio
```

This limits R's virtual memory to 16GB, adjust as needed, preventing system crashes during memory-intensive operations like parallel sequence alignment.

### Key Dependencies

**R Packages**:
- **Bioconductor**: `Biostrings`, `dada2`, `DECIPHER`, `pwalign`
- **Sequence Analysis**: `ape`, `seqinr`  
- **Data Analysis**: `tidyverse`, `dplyr`, `stringr`
- **Visualization**: `ggplot2`, `ComplexHeatmap`, `dendextend`
- **Parallel Processing**: `doSNOW`, `foreach`

**Python Packages** (for BLAST scripts):
- `biopython`, `pandas`, `pathlib`

**External Tools**:
- BLAST+ suite (`blastn`)
- VSEARCH (optional, for SINTAX classification)
