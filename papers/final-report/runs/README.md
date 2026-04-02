# Optimizer comparison runs

Each run of the full study should write a **self-contained folder** here (or anywhere you set `OPT_COMPARISON_OUTPUT_DIR`) with:

- `run_manifest.txt` — timestamp, git commit, Julia version, thread count, evaluation budget
- `comparison_results.txt` — text summary
- `optimizer_summary.csv` and `convergence_*.csv` — numeric traces
- `optimizer_comparison_table.tex` — LaTeX snippet for the report
- `optimizer_convergence_comparison.png`, `optimizer_best_so_far_convergence.png` — figures
- `best_optimized_configuration.mp4` — optional animation (large; often kept local only). Omit by setting `SKIP_COMPARISON_ANIMATION=1` before running the comparison.

## Preferred command

From the repository root:

```bash
chmod +x scripts/run_full_study.sh   # once
./scripts/run_full_study.sh --evals 50
```

Or with environment variables:

```bash
export JULIA_NUM_THREADS=8
MAX_EVALUATIONS=2 OPT_COMPARISON_OUTPUT_DIR=/tmp/opt_smoke julia --project=. -e 'include("scripts/optimization/compare_optimizers.jl")'
```

Large video files are ignored by `.gitignore` in this directory; commit small CSV/TeX/PNGs as needed for the report.
