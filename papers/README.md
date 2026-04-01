# Papers Directory

This directory contains LaTeX papers and research updates documenting the progress of this research project.

## Organization

The papers directory is organized as follows:

```
papers/
├── final-report/     # Final report (consolidated manuscript + PDF)
│   ├── Final_Report.tex
│   └── Final_Report.pdf
├── updates/          # Progress updates and interim reports
│   ├── Batra_Update_10_20_2025.tex
│   └── Batra_Update_7_23_25.pdf
├── technical/        # Standalone technical notes (e.g. spring force derivation)
│   └── spring_force_11x11.tex
├── figures/          # Figures and images used in papers
│   └── initial_state.png
├── build/            # Build artifacts (aux, log, etc.) - gitignored
└── README.md         # This file
```

## Naming Convention

Use descriptive names with dates:
- **`final-report/Final_Report.tex`** — primary final manuscript (PDF: `Final_Report.pdf`)
- `Batra_Update_MM_DD_YYYY.tex` for dated progress updates to Dr. Batra
- `update-YYYY-MM-DD.tex` for general progress updates
- `draft-YYYY-MM-DD.tex` for draft papers

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

**Important**: All build files (`.aux`, `.log`, `.out`, etc.) are automatically placed in the `build/` directory to keep source directories clean.

### Recommended: Use the compilation script

From the `papers/` directory:

```bash
./compile.sh final-report/Final_Report.tex
./compile.sh updates/Batra_Update_10_20_2025.tex
```

The script places the output PDF in the **same folder as the `.tex` file** (e.g. `final-report/` or `updates/`). Build artifacts go to `build/`.

Technical notes under `technical/` are usually compiled with `latexmk` directly (see below). From the repository root:

```bash
cd papers/technical && latexmk -pdf -output-directory=../build -aux-directory=../build spring_force_11x11.tex
```

The PDF is written to `papers/build/` unless you copy it next to the `.tex` source yourself.

### Manual compilation with latexmk

From the `papers/updates/` directory:

```bash
cd papers/updates
latexmk -pdf -output-directory=../build -aux-directory=../build Batra_Update_10_20_2025.tex
```

The `-output-directory` and `-aux-directory` flags ensure all build files go to `build/`.

### Manual compilation with pdflatex

From the `papers/updates/` directory:

```bash
cd papers/updates
pdflatex -output-directory=../build Batra_Update_10_20_2025.tex
pdflatex -output-directory=../build Batra_Update_10_20_2025.tex  # Run twice for references
```

### Using VS Code LaTeX Workshop (Build Button)

If you're using VS Code with the LaTeX Workshop extension, the build button (green play button) is configured to:
- Send all build artifacts (`.aux`, `.log`, etc.) to `papers/build/`
- Optionally copy the generated PDF next to the `.tex` file after successful compilation

The configuration is in `.vscode/settings.json` and includes:
- Output directory set to `../build` relative to the `.tex` file
- Auxiliary directory set to match output directory
- Post-build behavior for PDF placement (adjust if you build from `final-report/` vs `updates/`)

**Note**: The `.latexmkrc` files under `papers/`, `papers/updates/`, and `papers/final-report/` point outputs at `papers/build/`. If your editor copies PDFs only to `updates/`, update it when working in `final-report/`, or use `./compile.sh final-report/Final_Report.tex` from the command line.

### Compiling with figure paths:

When compiling from subdirectories, ensure figure paths are correct:
- From `updates/`: Use `../figures/figure_name.png`
- From root `papers/`: Use `figures/figure_name.png`

## Notes

- Keep source `.tex` files in version control
- Generated `.pdf` files are automatically placed in `updates/` directory and should be committed for easy access
- Build artifacts (`.aux`, `.log`, `.out`, etc.) are gitignored and stored in `build/`
- Use consistent citation styles (e.g., IEEE, APA, or journal-specific)
- For Unicode characters in verbatim environments, use ASCII equivalents (e.g., `m^-1` instead of `m⁻¹`)

## Common Issues

### Unicode Characters in Verbatim
If you encounter Unicode errors when compiling, replace Unicode characters in `verbatim` environments with ASCII equivalents:
- `m⁻¹` → `m^-1`
- `°` → `deg` or use `\textdegree`

### Figure Paths
Ensure figure paths are relative to the compilation directory. If compiling from `updates/`, use `../figures/`.
