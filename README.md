# pln_eval

Benchmark code and analysis for the paper:

**"When Does Poisson Log-Normal Modeling Help? A Systematic Comparison with Penalized Poisson Regression for Microbiome Count Prediction and Network Inference"**

## Overview

This repository provides the full pipeline for a LOTO-CV (Leave-One-Taxon-Out cross-validation) benchmark comparing PLN and GLMNet(Poisson) across 20 real microbiome count datasets and five network-inference benchmarks.

The main finding is that **N/D ratio** (samples-to-taxa) is the strongest predictor of which method wins, with **mean absolute correlation (MAC)** as the strongest secondary signal and overdispersion providing further discriminative power (combined AUC = 0.84).

## Paper

The manuscript is in `paper/archive/`.
Build with:

```bash
module load texlive/20250308
cd paper/archive
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Repository Structure

```
pln_eval/
├── bmr/                  # Benchmark result (.RData) files, one per dataset
├── out/                  # Computed outputs (CSV tables, metrics)
├── paper/
│   └── archive/          # LaTeX source for the arXiv submission
│       ├── main.tex
│       ├── supplement_content.tex
│       └── references.bib
├── plot_scripts/         # R scripts for generating figures
├── figures/              # Generated figures (PNG/PDF)
└── extract_se_pvalues.R  # Extracts SE and p-values from benchmark objects
```

## Datasets

20 real microbiome count datasets, sorted by N/D ratio (see Web Table 1 in the supplement).
The benchmark spans human gut, oral, soil, and other environments from publicly available studies.

## Methods Compared

**Count prediction:**
- `PLN` — Poisson Log-Normal (via the `PLNmodels` R package)
- `GLMNet(Poisson)` — penalized Poisson regression
- `Featureless` — baseline (mean predictor)

**Network inference:**
- `PLNNetwork` (sparse precision matrix via graphical lasso)
- `GLMNet(Poisson)` (nodewise regression)

## Key Results

- PLN wins when N/D is low and inter-taxon dependence (MAC) is high
- GLMNet(Poisson) is preferred as N/D grows
- For network inference, PLNNetwork is stronger on broad undirected community graphs; GLMNet(Poisson) is better matched to local or directional effects
- PLN is substantially more expensive computationally but remains practical at benchmark sizes

## Requirements

- R with `PLNmodels`, `glmnet`, `mlr3`, `atime` packages
- SLURM cluster (all R scripts should be run via `srun`)

```bash
srun /projects/genomic-ml/da2343/soak-r/bin/Rscript <script.R>
```
