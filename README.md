# Microbial Analysis Pipelines

This repository is associated with the work done in the IRTLab at the Laboratory of Microbial Genetics and Evolution, Graduate School of Life Sciences, Tohoku University.

This repository contains scripts, data structures, and workflows for analyzing microbial data generated from GC-MS, Sanger sequencing, next-generation sequencing (NGS), and growth curve experiments. The project aims to characterize microbial isolates, quantify phenanthrene degradation, and cluster microbial communities based on taxonomy and sequence similarity.

------------------------------------------------------------------------

## Project Structure

`project-root/├── data/ # Raw and intermediate input data`

`├── results/ # Output from scripts and analyses`

`├── scripts/ # All R, Python, and Quarto analysis scripts`

`├── databases/ # External databases (e.g. SILVA, EZBioCloud)`

`├── README.md # This file`

`├── .gitignore # Ignored files and folders`

`└── environment files # e.g. environment.yml or renv.lock`

------------------------------------------------------------------------

## Contents

### 1. Isolate Comparison (GC-MS)

**Goal**: Quantify phenanthrene degradation using GC-MS and visualize results.

-   `data/gcms_raw/`: `.txt` files with GC-MS outputs
-   `data/gcms_raw/`: `.xlsx` processed `.txt` to calculate phenanthrene/decane ratio and only save rows of interest
-   `scripts/gcms/plot_gc_data.R`: Loads processed `.xlsx`, plots phenanthrene/decane ratios or relative concentrations
-   Output: `results/gcms_plots/`

------------------------------------------------------------------------

### 2. Sanger Sequencing Pipeline

**Goal**: Process `.ab1` files, generate consensus sequences, identify taxa, and perform clustering.

-   `scripts/sanger/sanger_pipeline.py`: Trims and aligns reads using MAFFT
-   `scripts/sanger/run_blast.py`: BLASTs consensus sequences against EZBioCloud or SILVA
-   `scripts/sanger/16s_cluster_kmers.R`: Clusters consensus sequences using k-mer distance
-   `scripts/sanger/16s_cluster_alignment.R`: Clusters using PID from Needleman-Wunsch and Smith-Waterman alignments
-   Output: `results/sanger_alignments/`, `results/sanger_blast/`

------------------------------------------------------------------------

### 3. NGS Analysis (16S rRNA)

**Goal**: Process Illumina reads with DADA2, assign taxonomy, visualize ASVs.

-   `scripts/ngs/dada2_analysis.qmd`: DADA2 pipeline from raw reads to ASV table and relative abundance
-   `scripts/ngs/assign_taxa.py`: Assigns taxonomy using BLAST
-   `scripts/ngs/visualize_asvs.qmd`: Analyzes and visualizes ASV composition and clustering
-   `scripts/ngs/16s_asvs_vs_sanger.qmd`: visualizes the similarity between some ASVs and sanger alignments
-   Output: `results/ngs_taxa/`, `results/ngs_plots/`

> **Note**: Raw sequencing data is excluded from this repository. See `data/ngs/README.md` for how to organize your own data.

------------------------------------------------------------------------

### 4. Growth Curve Analysis

**Goal**: Fit growth models to OD600 data from a biophotorecorder.

-   `data/growth_curves/`: Contains `.dt3` and converted `.csv` files
-   `scripts/growth/growth_analysis.R`: Reads data, fits models, estimates `μ_max`, lag, generation time, and performs ANOVA
-   Output: `results/growth_models/`

------------------------------------------------------------------------

## Path Handling

All scripts assume you are running from the **repository root** and use **relative paths**.

------------------------------------------------------------------------

## Getting Started

### Setup Environment

Use R with `renv`, or Python with a virtual environment. To install dependencies:

### For Python

pip install -r requirements.txt

#### For R

renv::restore()

#### For Quarto

quarto install
