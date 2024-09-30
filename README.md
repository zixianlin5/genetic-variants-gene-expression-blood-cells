# Genetic Variants and Gene Expression in Blood Cells

## Overview

This repository hosts the source code and analytical scripts used in the project: **Investigating Genetic Variants Associated with Dynamic Gene Expression Profiles in Blood Cells**, which focuses on understanding the relationship between genetic variants and cell-state-dependent gene expression in blood cells using single-cell RNA sequencing data from the OneK1K cohort.

## Project Description

The goal of this study is to investigate the dynamic influence of genetic variants on gene expression in blood cells, incorporating cell-state-specific effects through expression principal components (ePCs). The analysis uses three Poisson mixed-effects (PME) regression models to assess the impact of genetic variants, and a trajectory inference framework is implemented to delineate the differentiation pathways. Additionally, the repository includes methods for identifying differential gene expression along inferred pseudotime trajectories.

### Key Components:

1. **Data Preprocessing**: Scripts for handling single-cell RNA-seq data and quality control.
2. **Regression Modeling**: Implementations of three PME models:
   - **Model 1**: Baseline gene expression model.
   - **Model 2**: Gene expression model with genetic variants (eSNPs).
   - **Model 3**: Model incorporating interactions between eSNPs and cell states (ePCs).
3. **Trajectory Analysis**: Methods for pseudotime inference using the Monocle3 package and iterative downsampling techniques, and procedures for combining results across multiple iterations to identify robust differentially expressed genes.

## Usage

1. **Data Preparation**: Use code `01 preprocess_onek1k.R` to `06 generate_genotype_matrix.R` to preprocess OneK1K data and prepare genotype matrix for the following steps.
2. **Model Fitting and Selection**: Use `07 fit_models.R` to fit models, `08 model_selection.R` to select model for each gene, and `09 generate_model_coeffs.R` to generate the selected models.
3. **Trajectory Analysis**: Use `10 convert_sparseMatrix.py` and `11 create_cds_obj.R` to prepare data, and `12 trajectory_analysis.R` to perform trajectory analysis on gene expression patterns over time using graph-based methods and conduct differential expression analysis.
4. **Visualization and Result Analysis**: The remaining unnumbered code is used for visualization and analysis of the results obtained from the previous code.
   
