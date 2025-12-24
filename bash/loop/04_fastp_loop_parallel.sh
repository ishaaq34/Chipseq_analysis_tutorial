#!/bin/bash
set -euo pipefail

# Requires parallel.py from fastp (often included or downloaded separately)
echo "Running fastp parallel.py script..."

# This assumes parallel.py is in the path or same directory
# python parallel.py -i /fastq_raw -o /fastq_cleaned -r /fastp_reports -f 3 -t 2
echo "Note: This requires 'parallel.py' from fastp instructions."
