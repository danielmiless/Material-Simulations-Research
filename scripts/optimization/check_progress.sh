#!/bin/bash
# Quick script to check optimization progress

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/comparison_output"

echo "============================================================"
echo "Optimization Progress Check"
echo "============================================================"
echo ""

# Check if processes are running
echo "1. Checking for running processes..."
JULIA_PROCS=$(ps aux | grep -E "julia.*optimization|run_overnight" | grep -v grep)
if [ -n "$JULIA_PROCS" ]; then
    echo "   ✓ Julia processes found:"
    echo "$JULIA_PROCS" | while read line; do
        echo "     $line"
    done
else
    echo "   ✗ No optimization processes currently running"
fi
echo ""

# Check screen sessions
echo "2. Checking screen sessions..."
SCREEN_SESSIONS=$(screen -ls 2>/dev/null | grep optimization)
if [ -n "$SCREEN_SESSIONS" ]; then
    echo "   ✓ Screen session found:"
    echo "$SCREEN_SESSIONS"
else
    echo "   ✗ No 'optimization' screen session found"
fi
echo ""

# Check for output directories
echo "3. Checking for output directories..."
if [ -d "$OUTPUT_DIR" ]; then
    OVERNIGHT_RUNS=$(ls -td "$OUTPUT_DIR"/overnight_run_* 2>/dev/null | head -1)
    if [ -n "$OVERNIGHT_RUNS" ]; then
        echo "   ✓ Latest overnight run: $(basename $OVERNIGHT_RUNS)"
        if [ -f "$OVERNIGHT_RUNS/overnight_run.log" ]; then
            echo "   Latest log entries:"
            tail -5 "$OVERNIGHT_RUNS/overnight_run.log" | sed 's/^/     /'
        fi
    else
        echo "   ✗ No overnight run directories found"
    fi
    
    # Check for regular comparison output files
    if [ -f "$OUTPUT_DIR/comparison_results.txt" ]; then
        echo "   ✓ Comparison results file found"
    fi
    if [ -f "$OUTPUT_DIR/optimizer_summary.csv" ]; then
        echo "   ✓ Optimizer summary CSV found"
    fi
else
    echo "   ✗ Output directory does not exist yet"
fi
echo ""

# Check log files
echo "4. Checking for log files..."
if [ -f "overnight_output.log" ]; then
    echo "   ✓ Main log file found (overnight_output.log)"
    echo "   Latest entries:"
    tail -5 overnight_output.log | sed 's/^/     /'
else
    echo "   ✗ No main log file found"
fi
echo ""

echo "============================================================"
echo "To reattach to screen session:"
echo "  TERM=xterm-256color screen -r optimization"
echo ""
echo "To view latest log:"
echo "  tail -f scripts/optimization/comparison_output/overnight_run_*/overnight_run.log"
echo "============================================================"

