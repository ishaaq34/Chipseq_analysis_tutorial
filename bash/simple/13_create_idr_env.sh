#!/bin/bash
set -euo pipefail

# Create IDR environment
if [ -f "idr_env.yml" ]; then
    conda env create -f idr_env.yml
else
    echo "Error: idr_env.yml not found."
    exit 1
fi

echo "To activate: conda activate idr_env"
