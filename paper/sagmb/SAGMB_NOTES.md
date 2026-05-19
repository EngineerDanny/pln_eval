# SAGMB Submission Package

This folder adapts the manuscript for Statistical Applications in Genetics and Molecular Biology using the De Gruyter LaTeX template downloaded from the official SAGMB page.

## Primary Files

- `main.tex`: SAGMB/De Gruyter-formatted manuscript.
- `main.pdf`: compiled manuscript PDF.
- `supplement.tex` and `supplement_content.tex`: supplementary material.
- `supplement.pdf`: compiled supplementary PDF.
- `title_page.tex`: separate title/declarations page, retained in case the submission system requests it separately.
- `title_page.pdf`: compiled title/declarations PDF.
- `references.bib`: bibliography copied from the current manuscript source.
- `figures/`: figure files referenced by the manuscript and supplement.

## Template Files

- `dgruyter.sty`, `dgruyter.ist`, `dgruyter.xdy`
- `journal-article.tex`, `journal-article.pdf`
- `Manual_dgruyter-template.pdf`, `FAQs_dgruyter-template.pdf`
- `SAGMB_LaTeX-Template-for-Authors.zip`

## Build Commands

Run these from `paper/sagmb`:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build/sagmb main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build/sagmb-supplement supplement.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build/sagmb-title title_page.tex
```

The current PDFs were compiled successfully with BibTeX and no unresolved citations or references.
