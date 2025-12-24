#!/bin/bash

# IDR Installation Script
# Creates 'idr_env' conda environment matching your existing setup

echo "========================================="
echo "Installing IDR in Conda Environment"
echo "========================================="
echo ""

# Initialize conda
eval "$(/opt/anaconda3/bin/conda shell.bash hook 2>/dev/null)" || \
eval "$(/opt/anaconda3/bin/conda shell.zsh hook 2>/dev/null)"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found"
    echo "Please run: conda init zsh"
    exit 1
fi

echo "✓ Found conda: $(conda --version)"
echo ""

# Check if environment already exists
if conda env list | grep -q "idr_env"; then
    echo "⚠️  Environment 'idr_env' already exists!"
    echo ""
    read -p "Do you want to remove and recreate it? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n idr_env -y
    else
        echo "Installation cancelled."
        exit 0
    fi
fi

echo "Creating conda environment 'idr_env'..."
conda create -n idr_env -y python=3.9 \
    -c conda-forge

echo ""
echo "Activating idr_env..."
source activate idr_env

echo ""
echo "Installing IDR and dependencies..."
pip install idr numpy scipy matplotlib

echo ""
echo "Testing installation..."
idr --version

echo ""
echo "========================================="
echo "✅ IDR Installed Successfully!"
echo "========================================="
echo ""
echo "Environment: idr_env"
echo "Location: $(conda info --envs | grep idr_env | awk '{print $NF}')"
echo ""
echo "To use IDR:"
echo "  1. Activate: conda activate idr_env"
echo "  2. Run: idr --version"
echo "  3. Deactivate: conda deactivate"
echo ""
echo "Example IDR command:"
echo "  idr --samples rep1.narrowPeak rep2.narrowPeak \\"
echo "      --input-file-type narrowPeak \\"
echo "      --rank p.value \\"
echo "      --output-file idr_output.txt \\"
echo "      --plot"
echo ""
echo "See Tutorial 13 for complete examples."
