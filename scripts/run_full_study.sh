#!/usr/bin/env bash
# Full optimizer comparison with timestamped artifacts for the final report.
#
# Usage:
#   ./scripts/run_full_study.sh [N]                 # N = MAX_EVALUATIONS (positional)
#   ./scripts/run_full_study.sh --evals 50
#   ./scripts/run_full_study.sh --quick             # small eval count for smoke tests
#   ./scripts/run_full_study.sh --output-dir PATH   # default: papers/final-report/runs/run_YYYYMMDD_HHMMSS/
#
# Environment:
#   OPT_COMPARISON_OUTPUT_DIR     — overrides --output-dir if set before calling (this script exports it)
#   MAX_EVALUATIONS               — used when no positional / --evals given (default 50)
#   JULIA_NUM_THREADS             — set before invoking for parallel evals (recommended on multi-core)
#   SKIP_COMPARISON_ANIMATION=1   — skip best-configuration mp4 (faster end-to-end)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

QUICK=false
EVALS=""
OUT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --quick)
      QUICK=true
      shift
      ;;
    --evals|--evaluations)
      EVALS="${2:?}"
      shift 2
      ;;
    --output|--output-dir)
      OUT="${2:?}"
      shift 2
      ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
    *)
      if [ -z "$EVALS" ] && [[ "$1" =~ ^[0-9]+$ ]]; then
        EVALS="$1"
        shift
      else
        echo "Unexpected argument: $1" >&2
        exit 1
      fi
      ;;
  esac
done

if [ "$QUICK" = true ] && [ -z "$EVALS" ]; then
  EVALS=8
fi
if [ -z "$EVALS" ]; then
  EVALS="${MAX_EVALUATIONS:-50}"
fi

if [ -z "$OUT" ]; then
  if [ -n "${OPT_COMPARISON_OUTPUT_DIR:-}" ]; then
    OUT="$OPT_COMPARISON_OUTPUT_DIR"
  else
    OUT="$ROOT/papers/final-report/runs/run_$(date +%Y%m%d_%H%M%S)"
  fi
fi

mkdir -p "$OUT"
export OPT_COMPARISON_OUTPUT_DIR="$OUT"
export MAX_EVALUATIONS="$EVALS"

echo "Repository: $ROOT"
echo "Output directory: $OUT"
echo "MAX_EVALUATIONS: $EVALS"
echo "JULIA_NUM_THREADS (current shell): ${JULIA_NUM_THREADS:-unset}"

exec julia --project=. -e 'include("scripts/optimization/compare_optimizers.jl")'
