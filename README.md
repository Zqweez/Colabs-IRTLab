# COLABS-IRTLab: Microbial Analysis Pipelines

This repository is associated with the work done in the IRTLab at the Laboratory of Microbial Genetics and Evolution, Graduate School of Life Sciences, Tohoku University.

This repository contains comprehensive scripts, data structures, and workflows for analyzing microbial data generated from GC-MS, qPCR gene expression analysis, Sanger sequencing, next-generation sequencing (NGS), and growth curve experiments. The project characterizes microbial isolates, quantifies phenanthrene degradation through both chemical analysis (GC-MS) and qPCR of degradation genes, performs taxonomic identification, and clusters microbial communities based on sequence similarity and phylogenetic relationships.

------------------------------------------------------------------------

## Project Structure

```
project-root/
├── data/                    # Raw and processed input data
│   ├── gcms_raw/           # GC-MS chromatogram data (.txt, .xlsx, .pdf)
│   ├── growth/             # Growth curve data (.dt files, .csv, .txt)
│   ├── phe_growth/         # Phenanthrene growth experiments
│   │   └── qPCR/           # qPCR data for gene expression analysis
│   │       ├── nidA/       # nidA gene expression data
│   │       └── phnAc/      # phnAc gene expression data
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
│   ├── qPCR/              # qPCR gene expression analysis (R)
│   └── sanger/            # Sanger sequencing analysis (R, Python)
├── install_packages.R      # Automated R package installation
├── setup_renv.R           # Reproducible environment setup
├── R_ENVIRONMENT_SETUP.md  # Detailed environment guide
├── README.md              # This file
└── 00 Colabs-repo.Rproj   # RStudio project file
```

------------------------------------------------------------------------

## Analysis Workflows

### 1. GC-MS and qPCR Analysis (Integrated Phenanthrene Degradation Study)

**Goal**: Quantify phenanthrene degradation using GC-MS chromatography and analyze gene abundance using qPCR. Also using qPCR of 16S rRNA to quantify the growth and using an external standard converting the Cq values to cells/mL. These analyses are conducted on the same sample sets to provide complementary evidence of degradation activity.

**Data**:
- `data/gcms_raw/L-tubes.xlsx`: Processed phenanthrene/decane ratios
- `data/gcms_raw/naphtol-calibration.xlsx`: 1-Naphthol calibration data
- `data/gcms_raw/*.txt`: Raw chromatogram outputs
- `data/gcms_raw/*.pdf`: Chromatogram visualizations
- `data/phe_growth/qPCR/*.xlsx`: qPCR results including Cq values, standard curves, and amplification data
- `data/phe_growth/qPCR/nidA/`: nidA gene expression data (ring-hydroxylating dioxygenase)
- `data/phe_growth/qPCR/phnAc/`: phnAc gene expression data (phenanthrene dioxygenase)

**Scripts**:
- `scripts/gcms/GCMS-phe-nap.R`: Analyzes phenanthrene degradation across samples and dates
- `scripts/gcms/GCMS-nap-calibration.R`: 1-Naphthol calibration curve analysis
- `scripts/qPCR/plot-qPCR.R`: Visualization of 16S rRNA qPCR results.
- `scripts/qPCR/analyze-qPCR.R`: Comprehensive qPCR data analysis and Cq calculations for 16S rRNA
- `scripts/qPCR/Calculate-standard.R`: Standard curve generation and efficiency calculations for 16S rRNA
- `scripts/qPCR/phe-genes-qPCR.R`: Gene expression analysis for phenanthrene degradation genes

**Features**:
- Integrated analysis of chemical degradation (GC-MS) and gene expression (qPCR)
- Standard curve validation and efficiency calculations
- Multi-gene expression profiling (nidA, phnAc)
- Statistical analysis of expression levels across isolates

**Output**: `results/gcms_plots/` - Chemical degradation analysis, qPCR gene expression plots, and integrated degradation assessments

---

### 2. Growth Curve Analysis (Phenanthrene-Supplemented Media)

**Goal**: Analyze microbial growth in phenanthrene-supplemented media using OD600 measurements from biophotorecorder.

**Note**: Basic LB media growth curves have been removed from the final analysis to focus on phenanthrene-specific growth responses.

**Data**:
- `data/phe_growth/*.dt*`: Growth with phenanthrene supplementation
- `data/phe_growth/*.csv`: Converted growth data with headers

**Scripts**:
- `scripts/growth/plot-growth.R`: Basic growth curve visualization and statistical modeling for LB curves
- `scripts/growth/plot-growth-phe.R`: Phenanthrene-supplemented growth analysis

**Features**:
- Automatic time correction for midnight rollovers
- Sliding window smoothing
- Focus on phenanthrene-utilizing isolates
- (Growth parameter estimation (μ_max, lag time, generation time) this never got implemented so the code is commented out)

**Output**: `results/growth_models/`, `results/growth_phe/` - Growth curves

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
- `scripts/sanger/run_blast.py`: Automated BLAST analysis against multiple databases
- `scripts/sanger/sanger_mafft_align.py`: Multiple sequence alignment using MAFFT

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
- `scripts/ngs/isolate-ASV-blast.py`: BLAST analysis comparing isolate sequences to ASVs
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

## Path Handling

All scripts assume you are running from the **repository root** and use **relative paths**. 

When working in RStudio, ensure your working directory is set to the project root:
```r
getwd()  # Should show the repository root
# If not, use: setwd("path/to/Colabs-IRTLab")
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
   git clone https://github.com/Zqweez/Colabs-IRTLab.git
   cd Colabs-IRTLab
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

**Python Packages** (for BLAST and alignment scripts):
- `biopython`, `pandas`, `pathlib`

**External Tools**:
- BLAST+ suite (`blastn`)
- MAFFT (for multiple sequence alignment)
- VSEARCH (optional, for SINTAX classification)
