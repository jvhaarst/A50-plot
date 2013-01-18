#!/bin/bash
rm Rplots* 2>/dev/null

R --no-save < AssemblyStats.R

if [[ -f "Rplots.pdf" ]]; then
    echo "converting pdf to png"
    bash pdf2png.sh
fi
