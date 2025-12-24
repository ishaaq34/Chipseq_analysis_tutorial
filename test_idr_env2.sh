#!/bin/bash

# Test script for creating and testing idr_env2

echo "========================================="
echo "Creating idr_env2 from YAML file"
echo "========================================="
echo ""

# Initialize conda
eval "$(/opt/anaconda3/bin/conda shell.bash hook 2>/dev/null)" || \
eval "$(/opt/anaconda3/bin/conda shell.zsh hook 2>/dev/null)"

# Create environment from YAML
echo "Creating environment from idr_env2.yml..."
conda env create -f idr_env2.yml

echo ""
echo "Activating idr_env2..."
conda activate idr_env2

echo ""
echo "Testing IDR installation..."
idr --version

echo ""
echo "Listing installed packages..."
conda list

echo ""
echo "========================================="
echo "✅ idr_env2 Created Successfully!"
echo "========================================="
echo ""
echo "To use this environment:"
echo "  conda activate idr_env2"
echo "  idr --version"
echo ""
echo "To remove if needed:"
echo "  conda env remove -n idr_env2"
