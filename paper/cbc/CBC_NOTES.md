# Computational Biology and Chemistry Submission Notes

Sources checked on May 3, 2026:

- Guide for authors: https://www.sciencedirect.com/journal/computational-biology-and-chemistry/publish/guide-for-authors
- Journal insights: https://www.sciencedirect.com/journal/computational-biology-and-chemistry/about/insights
- Elsevier LaTeX instructions: https://www.elsevier.com/en-gb/researcher/author/policies-and-guidelines/latex-instructions
- Elsevier/CTAN `elsarticle` bundle: https://ctan.org/texarchive/macros/latex/contrib/elsarticle

Relevant current requirements:

- Article type: Research article.
- Peer review: single anonymized review, so author names and affiliations remain in the manuscript.
- LaTeX: accepted; editable `.tex` source files are requested.
- References: author-year citations; this package uses `elsarticle-harv.bst`.
- Highlights: required as a separate editable file.
- Graphical abstract: required as a separate upload.
- Subscription route: no publication fee charged to authors; optional open access APC listed separately by Elsevier.

Package contents:

- `main.tex`: main Elsevier manuscript using `elsarticle`.
- `supplement.tex` and `supplement_content.tex`: standalone supplementary material.
- `title_page.tex`: separate title/declaration file in case Editorial Manager requests it.
- `highlights.tex`: separate highlights file.
- `references.bib`, `elsarticle.cls`, and Elsevier `.bst` files.
- `figures/`: figure files copied from `paper/archive`.
- `graphical_abstract/README.md`: graphical abstract upload note.
