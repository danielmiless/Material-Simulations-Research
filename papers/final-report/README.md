# Final Report

This folder holds the consolidated **final report** for the spring-mass sandwich-plate study (formerly `papers/updates/Batra_Update_01_27_2026.{tex,pdf}`).

| File | Purpose |
|------|---------|
| `Final_Report.tex` | LaTeX source |
| `Final_Report.pdf` | Built PDF (regenerate after edits) |

Figures are shared from `papers/figures/` (paths in the `.tex` file use `../figures/`).

## Compile

From `papers/`:

```bash
./compile.sh final-report/Final_Report.tex
```

The PDF is written next to `Final_Report.tex` in this directory; build artifacts go to `papers/build/`.
