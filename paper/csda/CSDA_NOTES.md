# CSDA LaTeX Preparation Notes

This folder adapts the archive manuscript to Elsevier's CAS LaTeX format for
Computational Statistics & Data Analysis.

Official sources checked:

- CSDA guide for authors: https://www.sciencedirect.com/journal/computational-statistics-and-data-analysis/publish/guide-for-authors
- Elsevier LaTeX instructions: https://www.elsevier.com/en-gb/researcher/author/policies-and-guidelines/latex-instructions
- Downloaded template ZIP: https://assets.ctfassets.net/o78em1y1w4i4/5uFmLZJTPDMAUjFnHRpjj8/6f19a979146eb93263763d87a894ab0d/els-cas-templates.zip

The source uses `cas-sc.cls`, the single-column Elsevier CAS class, because it is
the safer review-submission format. The original manuscript structure and text
were preserved from `../archive/main.tex`; only the journal wrapper, author
metadata, keyword block, bibliography style, and CAS compatibility setup were
changed.

Before final submission, add the real figure files with these base names:

- `figure1_conceptual_protocol`
- `figure3_count_prediction_case_studies`
- `figure4_network_f1_benchmarks`
- `figure5_omm12_network_comparison`
- `figure6_computational_scaling`
- `figureS1_host_fitness_network_comparison`

Compile on Monsoon with:

```sh
module load texlive/20250308
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

Compilation was tested with temporary placeholder PDFs for the missing figures.
The placeholders and generated build products were removed afterward.
