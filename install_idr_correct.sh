#!/bin/bash

# Correct IDR Installation Script
# This matches your working idr_env setup

echo "========================================="
echo "Installing IDR (ENCODE ChIP-seq Tool)"
echo "========================================="
echo ""

# Initialize conda
eval "$(/opt/anaconda3/bin/conda shell.bash hook 2>/dev/null)" || \
eval "$(/opt/anaconda3/bin/conda shell.zsh hook 2>/dev/null)"

# Create environment
echo "Creating conda environment 'idr_env'..."
conda create -n idr_env -y python=3.9 -c conda-forge

# Activate environment
echo "Activating idr_env..."
conda activate idr_env

# Install IDR from GitHub (correct source)
echo ""
echo "Installing IDR from source..."
pip install git+https://github.com/nboley/idr.git

# Or alternatively, install dependencies first then IDR
# pip install numpy scipy matplotlib
# pip install git+https://github.com/nboley/idr.git

# Test installation
echo ""
echo "Testing IDR installation..."
idr --version

echo ""
echo "========================================="
echo "✅ IDR Installed Successfully!"
echo "========================================="
echo ""
echo "IDR version should be 2.0.3 or later"
echo ""
echo "To use:"
echo "  conda activate idr_env"
echo "  idr --samples rep1.narrowPeak rep2.narrowPeak \\"
echo "      --input-file-type narrowPeak \\"
echo "      --rank p.value \\"
echo "      --output-file output.txt \\"
echo "      --plot"
