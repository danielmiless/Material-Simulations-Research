#!/bin/bash
# Script to copy PDF from build/ to updates/ directory
# This is called after LaTeX compilation

# Get the script directory (papers/)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="$SCRIPT_DIR/build"
UPDATES_DIR="$SCRIPT_DIR/updates"

# Get the PDF filename from the first argument (full path)
PDF_PATH="$1"

if [ -z "$PDF_PATH" ]; then
    echo "Error: No PDF path provided"
    exit 1
fi

# Extract just the filename
PDF_NAME=$(basename "$PDF_PATH")

# Copy PDF to updates directory
if [ -f "$PDF_PATH" ]; then
    mkdir -p "$UPDATES_DIR"
    cp "$PDF_PATH" "$UPDATES_DIR/"
    # Remove PDF from build directory after copying
    rm "$PDF_PATH"
    echo "Copied PDF to: $UPDATES_DIR/$PDF_NAME and removed from build directory"
else
    echo "Warning: PDF not found at $PDF_PATH"
    exit 1
fi

