# 5hmC Pathway Analysis Pipeline 🧬

A comprehensive toolkit for analyzing 5-hydroxymethylcytosine (5hmC) pathways in epigenetic modifications, featuring automated workflows for data processing, statistical analysis, and visualization.

## 📋 Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Pipeline Components](#pipeline-components)
- [Detailed Workflows](#detailed-workflows)
- [Output Structure](#output-structure)
- [External Resources](#external-resources)
- [Contributing](#contributing)

## Overview

This repository provides an end-to-end solution for analyzing 5-hydroxymethylcytosine (5hmC) data, with particular focus on:
- Raw sequencing data processing
- Differential expression analysis
- Feature counting and genomic annotations
- Machine learning with nested cross-validation
- Pathway enrichment analysis

# 5hmC Pathway Analysis Pipeline - Detailed Components

## 1. Raw Data Processing (`hmCSeal_v5.sh`)

### Inputs
- Raw FASTQ files (paired-end sequencing data)
- Sample names from `full_dataset_sample_names.txt` or `test_dataset_sample_names.txt`
- Reference genome files
- Annotation files for feature counting

### Processing Steps
1. **Trimming**: Uses trim_galore
   - Generates trimming reports
   - Outputs trimmed FASTQ files

2. **Alignment**: Uses bowtie2
   - Aligns to reference genome
   - Generates SAM files
   - Converts to BAM format

3. **BAM Processing**: Uses samtools
   - Sorts BAM files
   - Removes duplicates
   - Creates index files

4. **Feature Counting**: Uses featureCounts
   - Processes aligned BAM files
   - Counts reads based on annotations

### Outputs
- Trimmed FASTQ files
- Aligned and sorted BAM files (.bam)
- BAM index files (.bai)
- Three count matrices:
  - `ReadsCount_AutosomeGenes_ProteinCoding_v4_091823.txt`
  - `ReadsCount_AutosomeGenes_promoter_and_enhancer_v4_091823.txt`
  - `ReadsCount_AutosomeGenes_macs2_pe5_040124.txt`
- Log file: `hmCSeal_bowtie2_runv4_091823.log`

## 2. Differential Expression Analysis (`DESeq2_analysis.R`)

### Inputs
- Count matrices from feature counting step
- Metadata file with sample information
- Design formula for comparisons
- Cutoff parameters for significance

### Processing Steps
1. **Data Preprocessing**
   - Filtering low-count genes
   - Normalization
   - Sample quality control

2. **Statistical Analysis**
   - Differential expression testing
   - Multiple testing correction
   - Shrinkage estimation

### Outputs
- Normalized count matrices
- Differential expression results tables
- Visualizations in `/Output/DESeq2/`:
  - Heatmaps/
    - Expression heatmaps
    - Sample clustering
  - PCA/
    - PCA plots
    - Sample distribution
  - UMAP/
    - UMAP dimensionality reduction plots
- Statistical summary files
- Volcano plots
- Gene set enrichment results

## 3. Feature Analysis (`feature_count_fig.R`)

### Inputs
- BAM files from alignment
- Genomic feature annotations
- TxDb database
- Metadata file

### Processing Steps
1. **Feature Definition**
   - Promoters
   - UTRs (3' and 5')
   - Exons
   - Introns

2. **Count Generation**
   - Overlap counting
   - Feature distribution analysis

### Outputs
- Feature count matrices
- Distribution plots
- Statistical summaries
- Genomic region annotations
- Visualization plots:
  - Feature distribution boxplots
  - Count density plots

## 4. DHMR Analysis (`dhmr_macs_10_17_24.R`)

### Inputs
- Processed BAM files
- Peak calling parameters
- Sample metadata
- Genomic annotations

### Processing Steps
1. **Peak Calling**: Uses MACS3
2. **Consensus Peak Identification**: Uses DiffBind
3. **Statistical Analysis**
4. **Annotation**: Uses ChIPseeker

### Outputs
- Peak files
- Consensus region lists
- Statistical reports
- Annotated peaks
- Visualization plots

## 5. Machine Learning Analysis

### A. Nested Cross-validation (`nestedcv_glmnet_v3.R`)

#### Inputs
- Normalized count data
- Sample metadata
- Model parameters
- Cross-validation settings

#### Outputs
- Model performance metrics
- Feature importance scores
- ROC curves
- Cross-validation results
- Prediction probabilities

### B. Repeated Analysis (`repeated_nested_glmnet_v1.R`)

#### Inputs
- Results from nested cross-validation
- Multiple random seeds
- Validation parameters

#### Outputs
- Aggregated performance metrics
- Stability measures
- Comprehensive reports
- Best model selection results
- Feature importance consensus