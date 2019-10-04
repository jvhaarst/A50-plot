#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# This scripts calculates some genome statistitics, and generates a so called A50 plot.
# pylint: disable=line-too-long
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=E0401,E1101,C0111
# Needed modules
from __future__ import division, print_function
import os
import csv
import sys
import time
import itertools

# Handle compressed files
# From http://stackoverflow.com/a/26986344/194447
import gzip
import bz2

# Graphing stuff
import matplotlib
matplotlib.use('cairo') # Remove if using interactively (no outfile below)
matplotlib.use('Agg') # Remove if using interactively (no outfile below)
import matplotlib.pyplot as plt

def open_by_suffix(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rb')
    elif filename.endswith('.bz2'):
        return bz2.BZ2file(filename, 'r')
    else:
        return open(filename, 'r')

# Global variables
max_contigs = int(os.getenv('MAX_CONTIGS', 10000))  # Set to zero to ignore
min_length = int(os.getenv('MIN_LENGTH', 100))  # Set to zero to ignore
# Set the expected genome size in order to calculate the NG50, if zero, ignore
expected_genome_size = int(os.getenv('EXPECTED_GENOME_SIZE', 0)) # 0.350e9

DPI = 300
TITLE = os.getenv('TITLE', time.strftime("%Y%m%d"))
TYPE = os.getenv('TYPE', 'scaffold')

assemblies = {}
outfile = "assemblies_%i_%i_%s.png" % (min_length, max_contigs, TITLE)  # Leave empty if you want to use the interactive output instead of a file.
#outfile='' # Also outcomment/remove matplotlib.use above

# Colors to cycle through
# Adapted from http://colorbrewer2.org/?type=qualitative&scheme=Paired&n=12
# Usage from http://stackoverflow.com/questions/12236566/setting-different-color-for-each-series-in-scatter-plot-on-matplotlib
color_definition = [
    '#1f78b4',
    '#e31a1c',
    '#33a02c',
    '#ff7f00',
    '#a6cee3',
    '#fb9a99',
    '#b2df8a',
    '#fdbf6f',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928'
    ]

# http://stackoverflow.com/questions/4995733/how-to-create-a-spinning-command-line-cursor-using-python
spinner = itertools.cycle(['—', '\\', '|', '/'])

# http://stackoverflow.com/questions/9157314/python-write-data-into-csv-format-as-string-not-file
def csv2string(data):
    import csv
    from StringIO import StringIO
    si = StringIO()
    cw = csv.writer(si)
    cw.writerow(data)
    return si.getvalue().strip('\r\n')

# Class definitions
class Assembly:
    """Class to hold the assembly information, and perform the calculations"""
    @classmethod
    def count_gc(cls, seq):
        gc = sum(seq.count(x) for x in ['G', 'C', 'g', 'c', 'S', 's'])
        return gc

    @classmethod
    def count_n(cls, seq):
        n = sum(seq.count(x) for x in ['N', 'n'])
        return n

    def __init__(self, name, filename, min_length):
        # from __future__ import division
        from itertools import islice
        from Bio import SeqIO
        import sys
        self.name = name
        self.filename = filename
        self.min_length = min_length
        self.sizes = []
        self.gc = []
        self.n = []
        sys.stderr.write("Reading in %s\n" % self.filename)
        # Open handle to be able to use file
        with open_by_suffix(self.filename) as handle:
            # For each contig/scaffold, calculate length, and sort large to small
            # Also count the GC and N content per contig
            for rec in SeqIO.parse(handle, "fasta"):
                # Spinner
                sys.stderr.write(spinner.next())
                sys.stderr.flush()
                if len(rec) >= self.min_length:
                    self.sizes.append(len(rec))
                    self.gc.append(self.count_gc(rec.seq))
                    self.n.append(self.count_n(rec.seq))
                # Spinner
                sys.stderr.write('\b')
        sys.stderr.write("\n")
        self.sizes.sort(reverse=True)
        # Calculate all the statistics we can already do
        self.total_length = sum(self.sizes)
        self.count = len(self.sizes)
        self.sum = self.total_length
        try:
            self.max = max(self.sizes)
            self.min = min(self.sizes)
        except: # there is no contig smaller than maximum size
            raise
        self.average = sum(self.sizes) / len(self.sizes)
        self.median = self.sizes[int(len(self.sizes) / 2)]
        self.gc_content = sum(self.gc) * 100 /self.sum
        self.gc_sum = sum(self.gc)
        self.n_content = sum(self.n) * 100 /self.sum
        self.n_sum = sum(self.n)
        # Store the first number on the size list as the first of the incremental list
        self.incremental_sizes_cumulative = []
        self.incremental_sizes_cumulative.append(self.sizes[0])
        self.NG50 = 0
        self.LG50 = 0
        self.N50 = 0
        self.L50 = 0
        self.N90 = 0
        self.L90 = 0
        self.N95 = 0
        self.L95 = 0
        self.counter_over_1000 = 1 if (self.incremental_sizes_cumulative[-1] > 1000) else 0
        self.counter_over_10000 = 1 if (self.incremental_sizes_cumulative[-1] > 10000) else 0
        # Now iterate over the sizes list, and incrementally add it to the incremental list
        for index, size in enumerate(islice(self.sizes, 1, None)):
            self.incremental_sizes_cumulative.append(self.incremental_sizes_cumulative[-1] + size)
            # While we iterate over the array, we can gather NG50, N50, LG50 and the L50 when we reach half of the total assembly
            if expected_genome_size != 0:
                if (self.incremental_sizes_cumulative[-1] > expected_genome_size / 2) and (self.NG50 == 0):
                    self.NG50 = size
                    # We need to add one to the index because we started to loop at 1, not 0
                    self.LG50 = index + 1
            else:
                self.NG50 = 0
                self.LG50 = 0
            if (self.incremental_sizes_cumulative[-1] > self.total_length / 2) and (self.N50 == 0):
                self.N50 = size
                # We need to add one to the index because we started to loop at 1, not 0
                self.L50 = index + 1
            if (self.incremental_sizes_cumulative[-1] > self.total_length * 0.90) and (self.N90 == 0):
                self.N90 = size
                # We need to add one to the index because we started to loop at 1, not 0
                self.L90 = index + 1
            if (self.incremental_sizes_cumulative[-1] > self.total_length * 0.95) and (self.N95 == 0):
                self.N95 = size
                # We need to add one to the index because we started to loop at 1, not 0
                self.L95 = index + 1
            # Count the number of sequences larger than 1000
            if size >= 1000:
                self.counter_over_1000 += 1
            # Count the number of sequences larger than 1000
            if size >= 10000:
                self.counter_over_10000 += 1

    def return_stats(self):
        line = list()
        line.append(self.name)
        line.append(self.count)
        line.append(self.sum)
        line.append(self.max)
        line.append(self.min)
        line.append(self.average)
        line.append(self.median)
        line.append(self.N50)
        line.append(self.L50)
        line.append(self.NG50)
        line.append(self.LG50)
        line.append(self.N90)
        line.append(self.L90)
        line.append(self.N95)
        line.append(self.L95)
        line.append(self.counter_over_1000)
        line.append(self.counter_over_10000)
        line.append(self.gc_sum)
        line.append(self.gc_content)
        line.append(self.n_sum)
        line.append(self.n_content)
        return line

    def print_stats(self):
        # For prettyprinting the numbers
        import locale
        locale.setlocale(locale.LC_ALL, '')  # empty string for platform's default setting
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
        print("Name:%s") % self.name
        print("Count:%s") % format(self.count, "n")
        print("Sum:%s") % format(self.sum, "n")
        print("Max:%s") % format(self.max, "n")
        print("Min:%s") % format(self.min, "n")
        print("Average:%s") % format(self.average, "n")
        print("Median:%s") % format(self.median, "n")
        print("N50:%s") % format(self.N50, "n")
        print("L50:%s") % format(self.L50, "n")
        print("NG50:%s") % format(self.NG50, "n")
        print("LG50:%s") % format(self.LG50, "n")
        print("N90:%s") % format(self.N90, "n")
        print("L90:%s") % format(self.L90, "n")
        print("N95:%s") % format(self.N95, "n")
        print("L95:%s") % format(self.L95, "n")
        print("Count > 1000:%s") % format(self.counter_over_1000, "n")
        print("Count > 10000:%s") % format(self.counter_over_10000, "n")
        print("#GC:%s") % format(self.gc_sum, "n")
        print("GC percentage:%s") % format(self.gc_content, "n")
        print("#N:%s") % format(self.n_sum, "n")
        print("N percentage:%s") % format(self.n_content, "n")

# Check what kind of input we get
# Starts with ">" and is a single entry : 1 file
# Starts with ">" and there are more entries : multiple files
# Doesn't start with ">" : csv with 1 or more entries.
if len(sys.argv) == 1:
    sys.exit()
if len(sys.argv) == 2:
    # Single file, now check if the input is a FASTA formatted file
    with open_by_suffix(sys.argv[1]) as f:
        first_line = f.readline()
    if first_line[0] == '>':  # We got a single FASTA file
        input_file = sys.argv[1]
        assemblies[input_file] = Assembly(input_file, input_file, min_length)
    else:  # We assume a valid CSV file, put the values in the assemblies dict
        file_name = sys.argv[1]
        fp = open(file_name, 'rb')
        csv_file = csv.DictReader((row for row in fp if not row.startswith('#')), delimiter=',', quotechar='"')
        for row in csv_file:
            input_name = row['Seq_Name']
            input_file = row['Seq_File']
            assemblies[input_name] = Assembly(input_name, input_file, min_length)
if len(sys.argv) > 2:
    # Multiple files, for now assume they are all valid FASTA
    for input_file in sys.argv[1:]:
        assemblies[input_file] = Assembly(input_file, input_file, min_length)

# Now plot the A50 plots and print the stats
fig=plt.figure(dpi=DPI)
#A4
fig.set_size_inches(11.69,8.27)
fig.set_size_inches(8.27,11.69)
#A3
fig.set_size_inches(16.53,11.69)
fig.set_size_inches(11.69,16.53)

if TITLE != '': plt.suptitle(TITLE+"\nA50 plot of "+TYPE+" >"+str(min_length)+"bp", fontsize=20)

plt.subplot(2,1,1)
colors=itertools.cycle(color_definition)
print(csv2string(["Name", "Count", "Sum", "Max", "Min", "Average", "Median", "N50", "L50", "NG50", "LG50", "N90", "L90", "N95", "L95", "Count>1000", "Count>10000", "#GC", "GC", "#N", "N"]))
for name, assembly in iter(sorted(assemblies.items())):
    print(csv2string(assembly.return_stats()))
    color = next(colors)
    line = pylab.plot(assembly.incremental_sizes_cumulative, label=name, color=color)
    #line = pylab.plot(assembly.sizes, label=name+' sizes', color=color)
    #line = plt.plot(assembly.incremental_sizes, label=name, color=color)
plt.xlabel("Sequence count")
plt.ylabel("Incremental size (bp)")
plt.legend(loc='lower right')
if max_contigs > 0:
    plt.xlim([0, max_contigs])

plt.subplot(2,1,2)
colors=itertools.cycle(color_definition)
for name, assembly in iter(sorted(assemblies.items())):
    color = next(colors)
    line = plt.plot(assembly.incremental_sizes,assembly.sizes, label=name, color=color)
plt.xlabel("Incremental size (bp)")
plt.ylabel(TYPE+" size (bp)")
plt.legend(loc='upper right')
if max_contigs > 0:
    plt.xlim([0, assembly.incremental_sizes[max_contigs]])

plt.tight_layout()
fig.subplots_adjust(top=0.8)

if outfile != '':
    sys.stderr.write("Writing %s\n" % outfile)
    plt.savefig(outfile, dpi=(DPI))
else:
    plt.show()
