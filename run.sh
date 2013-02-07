#!/bin/bash
INFILE=input.csv
rm $INFILE.* 2>/dev/null

Rscript AssemblyStats.R $INFILE

