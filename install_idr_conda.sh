#!/bin/bash

# IDR Installation Script - Conda Environment Version
# Creates a dedicated conda environment for IDR

echo "========================================="
echo "Installing IDR in Conda Environment"
echo "========================================="
echo ""

# Initialize conda for this session
eval "$(conda shell.bash hook 2>/dev/null)" || eval "$(conda shell.zsh hook 2>/dev/null)" || true

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found"
    echo "Please initialize conda first:"
    echo "  conda init zsh  # or bash"
    echo "  source ~/.zshrc # or ~/.bashrc"
    exit 1
fi

echo "Found conda: $(which conda)"
echo "Conda version: $(conda --version)"
echo ""

# Create idr_env environment
echo "Creating conda environment 'idr_env'..."
conda create -n idr_env -y python=3.9

# Activate the environment
echo ""
echo "Activating idr_env..."
conda activate idr_env

# Install IDR
echo ""
echo "Installing IDR..."
pip install idr

# Test installation
echo ""
echo "Testing IDR installation..."
idr --version

echo ""
echo "✅ IDR installed successfully in 'idr_env' environment!"
echo ""
echo "========================================="
echo "How to Use IDR:"
echo "========================================="
echo ""
echo "1. Activate the environment:"
echo "   conda activate idr_env"
echo ""
echo "2. Run IDR:"
echo "   idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak \\"
echo "       --input-file-type narrowPeak \\"
echo "       --rank p.value \\"
echo "       --output-file idr_output.txt \\"
echo "       --plot"
echo ""
echo "3. Deactivate when done:"
echo "   conda deactivate"
echo ""
echo "See Tutorial 13 for complete examples."
