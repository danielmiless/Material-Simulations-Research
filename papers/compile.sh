#!/bin/bash
# Compilation script for LaTeX papers
# Ensures all build files go to the build/ directory

# Get the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="$SCRIPT_DIR/build"

# Create build directory if it doesn't exist
mkdir -p "$BUILD_DIR"

# Check if a .tex file was provided
if [ -z "$1" ]; then
    echo "Usage: ./compile.sh <paper-name.tex>"
    echo "Example: ./compile.sh updates/Batra_Update_10_20_2025.tex"
    exit 1
fi

# Get the paper name without extension
PAPER_PATH="$1"
PAPER_DIR=$(dirname "$PAPER_PATH")
PAPER_NAME=$(basename "$PAPER_PATH" .tex)

# Change to the directory containing the .tex file
cd "$SCRIPT_DIR/$PAPER_DIR" || exit 1

# Compile with output directory set to build/
echo "Compiling $PAPER_NAME.tex..."
echo "Build files will be placed in: $BUILD_DIR"

# Use latexmk with output directory
latexmk -pdf \
    -output-directory="$BUILD_DIR" \
    -aux-directory="$BUILD_DIR" \
    "$PAPER_NAME.tex"

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "Compilation successful!"
    echo "PDF location: $BUILD_DIR/$PAPER_NAME.pdf"
    echo ""
    echo "To copy PDF to updates directory:"
    echo "  cp $BUILD_DIR/$PAPER_NAME.pdf $SCRIPT_DIR/updates/"
else
    echo ""
    echo "Compilation failed. Check the log file: $BUILD_DIR/$PAPER_NAME.log"
    exit 1
fi

