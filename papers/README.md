# Papers Directory

This directory contains LaTeX papers and research updates documenting the progress of this research project.

## Organization

Please organize papers using the following structure:

```
papers/
├── updates/          # Progress updates and interim reports
├── drafts/          # Draft papers in preparation
└── final/            # Finalized papers (if any)
```

## Naming Convention

Use descriptive names with dates:
- `update-YYYY-MM-DD.tex` for progress updates
- `draft-YYYY-MM-DD.tex` for draft papers
- `paper-title.tex` for final papers

## LaTeX Template

When creating new papers, consider including:
- Title and author information
- Abstract
- Introduction
- Methodology
- Results
- Discussion
- References

## Compilation

Papers can be compiled using standard LaTeX tools:
```bash
pdflatex paper-name.tex
bibtex paper-name
pdflatex paper-name.tex
pdflatex paper-name.tex
```

Or use `latexmk` for automatic compilation:
```bash
latexmk -pdf paper-name.tex
```

## Notes

- Keep source `.tex` files in version control
- Generated `.pdf` files should be committed for easy access
- Use consistent citation styles (e.g., IEEE, APA, or journal-specific)

