# MHE_rif_R

Shotgun metagenomic analysis pipeline and downstream statistical workflows for patients with minimal hepatic encephalopathy (MHE) treated with rifaximin, including cohorts from the United Kingdom and Spain.

This repository contains scripts for:

· Raw metagenomic preprocessing and quality control
· Taxonomic and functional microbiome profiling
· Antibiotic resistance gene (ARG) analysis
· Metagenome-assembled genome (MAG) reconstruction and annotation
· Comparative microbiome statistics and visualization
· Longitudinal before/after rifaximin treatment analyses

The project integrates shotgun sequencing datasets from independent UK and Spanish patient cohorts to investigate microbiome, resistome, and strain-level changes associated with rifaximin therapy.

### 📁 Repository Structure
```
MHE_rif_R/
├── ARGs/
├── MAGs/
├── Microbiome_analysis_and_Statistics/
└── Pre_processing/
```

### 📘 Overview

#### 1. Pre-processing

This module includes:

· Quality control of raw sequencing reads
· Functional and taxonomic annotation using NGLess
· MetaPhlAn profiling
· StrainPhlAn marker extraction and phylogenetic analyses
· Environment configuration files

Requirements
NGLess
MetaPhlAn
StrainPhlAn
SLURM environment (recommended)

#### 2. ARG Analysis

This workflow identifies and quantifies antibiotic resistance genes from shotgun metagenomic assemblies and mappings.

· Assembly using MEGAHIT
· ARG annotation using:
    · ResFinder
    · RGI / CARD
· Read mapping to ARG-containing contigs
· ARG abundance normalization


#### 3. MAG Reconstruction and Annotation

This module reconstructs and characterizes metagenome-assembled genomes (MAGs).

· Contig binning
· MAG quality control
· Taxonomic classification
· Functional annotation
· Dereplication
· Phylogenomics
· KEGG/GMM functional analyses
· Pan-genome analyses


#### 4. Microbiome Analysis and Statistics

This section contains statistical analyses and figure-generation scripts.

· Phyloseq object construction
· Alpha diversity
· Beta diversity
· Differential abundance testing
· Batch correction
· Network analyses
· Circos plots
· Sankey diagrams
· MaAsLin2 modeling
· Longitudinal before/after treatment comparisons


The workflows were developed primarily for execution in a Linux HPC environment with SLURM scheduling.
