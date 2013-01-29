#!/bin/bash
rm Rplots* 2>/dev/null

Rscript AssemblyStats.R input.csv

if [[ -f "Rplots.pdf" ]]; then
    echo "converting pdf to png"
    bash pdf2png.sh
fi
