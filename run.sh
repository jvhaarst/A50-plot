#!/bin/bash
INFILE=$1
rm $INFILE.* 2>/dev/null

Rscript AssemblyStats.R $INFILE

#for FILE in *.csv; do rm $FILE.png; Rscript AssemblyStats.R $FILE &\ ; done
