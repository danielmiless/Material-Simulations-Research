#!/bin/bash
# Overnight Optimizer Comparison Script
# Runs the full comparison with higher max_evaluations and logs everything

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/comparison_output"
MAX_EVALS="${1:-50}"  # Default to 50 for overnight run, can override

# Create timestamped output directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="$OUTPUT_DIR/overnight_run_$TIMESTAMP"
mkdir -p "$RUN_DIR"

# Log file
LOG_FILE="$RUN_DIR/overnight_run.log"

echo "============================================================" | tee -a "$LOG_FILE"
echo "Overnight Optimizer Comparison Run" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "Start Time: $(date)" | tee -a "$LOG_FILE"
echo "Max Evaluations: $MAX_EVALS" | tee -a "$LOG_FILE"
echo "Output Directory: $RUN_DIR" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Function to run with logging (macOS-compatible, no timeout)
# Uses tee to show output in real-time AND log to file
run_with_logging() {
    local description=$1
    local command=$2
    
    echo "[$(date +'%H:%M:%S')] Starting: $description" | tee -a "$LOG_FILE"
    echo "Command: $command" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Run command and show output in real-time while also logging to file
    # Use tee to duplicate output to both stdout and log file
    bash -c "$command" 2>&1 | tee -a "$LOG_FILE"
    
    local exit_code=${PIPESTATUS[0]}  # Get exit code of the command, not tee
    if [ $exit_code -eq 0 ]; then
        echo "" | tee -a "$LOG_FILE"
        echo "[$(date +'%H:%M:%S')] Completed: $description" | tee -a "$LOG_FILE"
    else
        echo "" | tee -a "$LOG_FILE"
        echo "[$(date +'%H:%M:%S')] ERROR: $description (exit code: $exit_code)" | tee -a "$LOG_FILE"
    fi
    echo "" | tee -a "$LOG_FILE"
    return $exit_code
}

# Change to project root
cd "$PROJECT_ROOT"

# Run the comparison script with logging
# All optimizers run in parallel by default
run_with_logging "Full Optimizer Comparison" \
    "julia --project=. -e 'MAX_EVALUATIONS=$MAX_EVALS; include(\"$SCRIPT_DIR/compare_optimizers.jl\")'"

# Copy results to timestamped directory
# The comparison script saves to OUTPUT_DIR directly, not OUTPUT_DIR/comparison_output
if [ -f "$OUTPUT_DIR/comparison_results.txt" ]; then
    echo "Copying comparison results to run directory..." | tee -a "$LOG_FILE"
    cp "$OUTPUT_DIR"/comparison_results.txt "$RUN_DIR/" 2>/dev/null || true
    cp "$OUTPUT_DIR"/optimizer_summary.csv "$RUN_DIR/" 2>/dev/null || true
    cp "$OUTPUT_DIR"/optimizer_comparison_table.tex "$RUN_DIR/" 2>/dev/null || true
    cp "$OUTPUT_DIR"/optimizer_convergence_comparison.png "$RUN_DIR/" 2>/dev/null || true
    cp "$OUTPUT_DIR"/convergence_*.csv "$RUN_DIR/" 2>/dev/null || true
    cp "$OUTPUT_DIR"/best_optimized_configuration.mp4 "$RUN_DIR/" 2>/dev/null || true
fi

# Final summary
echo "" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "Overnight Run Complete" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "End Time: $(date)" | tee -a "$LOG_FILE"
echo "Results saved in: $RUN_DIR" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Optional: Send notification (uncomment and configure if desired)
# echo "Overnight optimization run completed at $(date)" | mail -s "Optimization Complete" your-email@example.com

