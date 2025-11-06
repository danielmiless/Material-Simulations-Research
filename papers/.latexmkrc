# LaTeXmk configuration for papers directory
# This ensures all build files go to the build/ directory
# Works for command-line compilation and IDE build buttons

# Set output directory for all build artifacts
# When compiling from updates/, this goes to ../build/
# When compiling from papers/, this goes to ./build/
$out_dir = '../build';

# Ensure PDFs are also placed in build directory
$pdf_mode = 5;  # Use pdflatex

# Set aux directory to match output directory
$aux_dir = $out_dir;

# Clean up build files
$clean_ext = 'aux bbl blg fdb_latexmk fls log out synctex.gz toc lof lot idx ilg ind nav snm vrb';

# Don't clean PDFs (we want to keep them)
$clean_full_ext = 'aux bbl blg fdb_latexmk fls log out synctex.gz toc lof lot idx ilg ind nav snm vrb';

