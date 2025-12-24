#!/bin/bash
set -euo pipefail

# Check conda version
echo "Checking Conda version..."
conda --version

# Create environment (assumes chip_env.yml exists)
if [ -f "chip_env.yml" ]; then
    echo "Creating environment from chip_env.yml..."
    conda env create -f chip_env.yml
else
    echo "Error: chip_env.yml not found. Please create it first."
    exit 1
fi

# Note: conda activate often requires 'source' in scripts or shell init
echo "To activate the environment, run: conda activate chip"

# Verification commands (run these after activating)
# which fastqc
# which bowtie2
# which macs3
# fastqc --version
# bowtie2 --version
