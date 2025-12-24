#!/bin/bash
set -euo pipefail

# Generate MACS3 model plots from all model.r files

echo "Generating MACS3 model plots..."

# Find and process all model.r files
for MODEL_FILE in macs3_results/*_model.r; do
    if [ -f "$MODEL_FILE" ]; then
        echo "Processing: $MODEL_FILE"
        
        # Run Rscript
        Rscript "$MODEL_FILE"
        
        # Get PDF name
        PDF_FILE="${MODEL_FILE%.r}.pdf"
        
        if [ -f "$PDF_FILE" ]; then
            echo "✅ Generated: $PDF_FILE"
        fi
    fi
done

echo "Done!"
