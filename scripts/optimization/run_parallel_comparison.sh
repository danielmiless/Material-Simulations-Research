#!/bin/bash
# Parallel Optimizer Comparison Script
# Runs each optimizer in a separate Julia process for parallel execution

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/comparison_output"
MAX_EVALS="${1:-20}"  # Default to 20, can override with first argument

echo "============================================================"
echo "Parallel Optimizer Comparison"
echo "============================================================"
echo "Max Evaluations: $MAX_EVALS"
echo "Output Directory: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create a temporary directory for individual results
TEMP_DIR="$OUTPUT_DIR/temp_parallel_$$"
mkdir -p "$TEMP_DIR"

echo "Starting optimizers in parallel..."
echo "Results will be saved to: $OUTPUT_DIR"
echo ""

# Function to run an optimizer
run_optimizer() {
    local optimizer_script=$1
    local optimizer_name=$2
    local log_file="$TEMP_DIR/${optimizer_name}.log"
    
    echo "[$(date +'%H:%M:%S')] Starting $optimizer_name..."
    
    cd "$PROJECT_ROOT"
    julia --project=. "$SCRIPT_DIR/$optimizer_script" "$MAX_EVALS" > "$log_file" 2>&1 &
    local pid=$!
    echo "$pid" > "$TEMP_DIR/${optimizer_name}_pid.txt"
    echo "[$(date +'%H:%M:%S')] $optimizer_name started (PID: $pid)"
}

# Run optimizers in parallel (background processes)
# Note: Each optimizer script should accept max_evaluations as first argument
run_optimizer "optimize_optim.jl" "optim" || true
run_optimizer "optimize_evolutionary.jl" "evolutionary" || true

# Metaheuristics - we'll run these separately since they need algorithm parameter
# For now, just run PSO as an example
echo "[$(date +'%H:%M:%S')] Note: Metaheuristics algorithms should be run via compare_optimizers.jl"
echo ""

echo "Waiting for optimizers to complete..."
echo "Monitor progress with: tail -f $TEMP_DIR/*.log"
echo ""

# Wait for all background processes
wait

echo ""
echo "============================================================"
echo "Parallel execution complete!"
echo "============================================================"
echo "Check logs in: $TEMP_DIR"
echo "For full comparison with all optimizers, use: julia compare_optimizers.jl $MAX_EVALS"
