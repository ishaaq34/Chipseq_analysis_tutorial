#!/bin/bash
set -euo pipefail

samtools flagstat ceb_ENCFF327JFG.bam > ceb_ENCFF327JFG.flagstat.txt
cat ceb_ENCFF327JFG.flagstat.txt
