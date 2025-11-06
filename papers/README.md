# Papers Directory

This directory contains LaTeX papers and research updates documenting the progress of this research project.

## Organization

The papers directory is organized as follows:

```
papers/
├── updates/          # Progress updates and interim reports
│   ├── Batra_Update_10_20_2025.tex
│   └── Batra_Update_7_23_25.pdf
├── figures/          # Figures and images used in papers
│   └── initial_state.png
├── build/            # Build artifacts (aux, log, etc.) - gitignored
└── README.md         # This file
```

## Naming Convention

Use descriptive names with dates:
- `Batra_Update_MM_DD_YYYY.tex` for progress updates to Dr. Batra
- `update-YYYY-MM-DD.tex` for general progress updates
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

**Important**: All build files (`.aux`, `.log`, `.out`, etc.) are automatically placed in the `build/` directory to keep source directories clean.

### Recommended: Use the compilation script

From the `papers/` directory:

```bash
./compile.sh updates/Batra_Update_10_20_2025.tex
```

This ensures all build artifacts go to `build/` automatically.

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

If you're using VS Code with the LaTeX Workshop extension, the build button (green play button) is configured to automatically send all build files to `papers/build/`. 

The configuration is in `.vscode/settings.json` and includes:
- Output directory set to `../build` relative to the `.tex` file
- Auxiliary directory set to match output directory
- Automatic cleanup of build files from source directories

**Note**: The `.latexmkrc` files in `papers/` and `papers/updates/` directories are configured to automatically use the `build/` directory. The VS Code settings work together with these to ensure build files always go to the correct location.

### Compiling with figure paths:

When compiling from subdirectories, ensure figure paths are correct:
- From `updates/`: Use `../figures/figure_name.png`
- From root `papers/`: Use `figures/figure_name.png`

## Notes

- Keep source `.tex` files in version control
- Generated `.pdf` files should be committed for easy access
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
