#!/bin/bash
rm Rplots.pdf  2>/dev/null
rm Rplots*.png 2>/dev/null

R --no-save < AssemblyStats.R

if [[ -f "Rplots.pdf" ]]; then
    bash pdf2png.sh
fi
