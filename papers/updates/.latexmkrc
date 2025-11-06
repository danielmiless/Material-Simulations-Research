# LaTeXmk configuration for updates subdirectory
# This ensures all build files go to the build/ directory

# Set output directory for all build artifacts
$out_dir = '../build';

# Set aux directory to match output directory  
$aux_dir = $out_dir;

# Ensure PDFs are also placed in build directory
$pdf_mode = 5;  # Use pdflatex

# Clean up build files
$clean_ext = 'aux bbl blg fdb_latexmk fls log out synctex.gz toc lof lot idx ilg ind nav snm vrb';

# Don't clean PDFs (we want to keep them)
$clean_full_ext = 'aux bbl blg fdb_latexmk fls log out synctex.gz toc lof lot idx ilg ind nav snm vrb';

