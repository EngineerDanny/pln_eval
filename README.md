# PLN Evaluation Figure Plan

This repository contains the benchmark and graph-comparison work for the PLN evaluation project.

The current paper figure plan is intentionally narrow. The figure set should match the story we can defend with the current results:

- broad predictive benchmarking across many microbiome datasets
- when PLN helps
- a small case-study benchmark figure with representative PLN wins
- comparison of the retained graph methods
- computational cost
- one representative network visualization

## Figures To Keep

### 1. Cross-Dataset Benchmark Overview

Purpose:
- show predictive performance across the completed datasets
- establish the main benchmark result across the broad dataset collection

Suggested format:
- heatmap, dot plot, or ranked interval plot
- one point or cell per dataset-method combination

Core inputs:
- [analysis/full_comparison_table.csv](/projects/genomic-ml/da2343/PLN/pln_eval/analysis/full_comparison_table.csv)
- dataset summaries in [bmr](/projects/genomic-ml/da2343/PLN/pln_eval/bmr)

Main methods:
- `Featureless`
- `GLMNet (Poisson)`
- `PLN`

Main message:
- across a diverse set of datasets, how often and by how much PLN improves predictive performance

### 2. Case-Study Benchmark Figure

Purpose:
- show a small set of representative datasets where PLN clearly beats the baselines
- keep the visual narrative simple and interpretable

Suggested format:
- one-row faceted benchmark figure
- x-axis: mean Poisson deviance
- methods shown:
  - `PLN`
  - `GLMNet (Poisson)`
  - `Featureless`
- annotations:
  - `PLN - GLMNet` difference
  - paired p-value from the saved benchmark object

Current script:
- [plot_case_study_benchmark.R](/projects/genomic-ml/da2343/PLN/pln_eval/plot_scripts/plot_case_study_benchmark.R)

Current figure:
- [pln_case_study_benchmark_4datasets.png](/projects/genomic-ml/da2343/PLN/pln_eval/figures/pln_case_study_benchmark_4datasets.png)

Current datasets in that figure:
- `amgut2`
- `crc_zeller`
- `hiv_lozupone_family`
- `cdi_schubert_family`

Main message:
- PLN can beat both `GLMNet (Poisson)` and `Featureless` on real benchmark datasets without relying on obviously extreme outlier panels

### 3. Dataset Characteristics Versus Method Advantage

Purpose:
- explain when PLN helps
- connect dataset structure to method performance

Suggested format:
- scatterplots
- optional faceting by taxonomic level or cohort family

Candidate x-axes:
- `p / n`
- zero fraction / sparsity
- covariance-related summaries such as MAC
- mean abundance variation

Candidate y-axis:
- PLN minus lasso improvement
- PLN minus featureless improvement

Core inputs:
- [analysis/covariance_metrics.csv](/projects/genomic-ml/da2343/PLN/pln_eval/analysis/covariance_metrics.csv)
- [analysis/full_comparison_table.csv](/projects/genomic-ml/da2343/PLN/pln_eval/analysis/full_comparison_table.csv)

Working interpretation from the current case-study datasets:
- higher MAC can favor PLN
- low `n / p` can favor PLN
- high sparsity can also favor PLN even when `n > p`
- so the useful rule is not just `p > n`; it is broader than that

Main message:
- characterize the regimes where PLN is most useful without reducing the result to one oversimplified condition

### 4. Graph-Method Comparison

Purpose:
- compare the retained graph methods under one common graph score
- support the decision to keep the LOTO-style PLN graph construction

Current retained methods from [network_loglik_test.R](/projects/genomic-ml/da2343/PLN/pln_eval/network_consistency/network_loglik_test.R):
- `Baseline_diagonal`
- `PLN_glasso`
- `LOTO_PLN_glasso`
- `LOTO_glmnet_CV1se`

Current score:
- held-out Gaussian pseudo-loglikelihood on `log1p` data

Suggested format:
- bar plot with replicate points
- or interval plot showing mean and replicate-level values

Main message:
- `LOTO_PLN_glasso` is currently the strongest graph method under the shared graph score

### 5. Runtime / Memory Figure

Purpose:
- show computational cost, not just accuracy
- make the tradeoff between graph quality and compute explicit

Suggested format:
- two-panel figure
- left: runtime
- right: memory

Methods to include:
- `PLN_glasso`
- `LOTO_PLN_glasso`
- `LOTO_glmnet_CV1se`

Possible sources:
- timings printed by [network_loglik_test.R](/projects/genomic-ml/da2343/PLN/pln_eval/network_consistency/network_loglik_test.R)
- SLURM accounting for memory if available

Main message:
- `LOTO_PLN_glasso` appears strongest, but it is also more expensive

### 6. Representative Network Figure

Purpose:
- include at least one interpretable network visualization
- show qualitative differences between the retained graph methods

Suggested dataset:
- `amgut2`

Suggested methods:
- `PLN_glasso`
- `LOTO_PLN_glasso`
- `LOTO_glmnet_CV1se`

Suggested format:
- 3 side-by-side network panels
- same node layout across methods
- plot only a capped edge set, such as top edges by absolute weight

Recommended visual encoding:
- edge color: sign
- edge width: absolute weight
- node size: abundance or degree

Main message:
- the different graph constructions do not produce the same network, and the LOTO-style PLN result should be visually interpretable

## Figures To Avoid For Now

These are intentionally not part of the current paper story:

- network consistency or stability figures from the deleted workflow
- standalone win-rate or paired-improvement plots
- sample-size scaling as a main-text figure
- exploratory diagnostics based on arbitrary fixed-parameter covariance constructions

## Current Method Scope

For predictive benchmarking:
- `Featureless`
- `GLMNet (Poisson)`
- `PLN`

For graph comparisons:
- `Baseline_diagonal`
- `PLN_glasso`
- `LOTO_PLN_glasso`
- `LOTO_glmnet_CV1se`

Files still relevant for graph work:
- [network_loglik_test.R](/projects/genomic-ml/da2343/PLN/pln_eval/network_consistency/network_loglik_test.R)

## Immediate Next Plotting Priorities

If plotting work resumes, the order should be:

1. cross-dataset benchmark overview
2. case-study benchmark figure refinement
3. dataset-characteristics versus method advantage
4. graph-method comparison
5. runtime/memory figure
6. representative network figure
