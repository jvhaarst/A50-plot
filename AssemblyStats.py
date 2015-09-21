#!/usr/bin/env python2
# This scripts calculates some genome statistitics, and generates a so called A50 plot.
# Needed modules
from __future__ import division
import csv
import sys
import pylab
import itertools

# Global variables
max_contigs = 0000  # Set to zero to ignore
min_length = 200  # Set to zero to ignore
# Set the expected genome size in order to calculate the NG50, if zero, ignore
expected_genome_size = 0 # 0.350e9

DPI = 300
TITLE = ''
assemblies = {}
outfile = "assemblies_%i_%s.png" % (max_contigs,TITLE)  # Leave empty if you want to use the interactive output instead of a file.
#outfile=''

# Colors to cycle through
# Adapted from http://colorbrewer2.org/?type=qualitative&scheme=Paired&n=12
# Usage from http://stackoverflow.com/questions/12236566/setting-different-color-for-each-series-in-scatter-plot-on-matplotlib
colors = itertools.cycle([
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
    ])

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
    def __init__(self, name, filename,min_length):
        # from __future__ import division
        from itertools import islice
        from Bio import SeqIO
        import sys
        self.name = name
        self.filename = filename
        self.min_length=min_length
        self.sizes = []
        sys.stderr.write("Reading in %s\n" % self.filename)
        # For each contig/scaffold, calculate length, and sort large to small
        for rec in SeqIO.parse(open(self.filename), "fasta"):
            if (len(rec) >= self.min_length):
                self.sizes.append(len(rec))
            else:
                if (self.min_length == 0):
                    self.sizes.append(len(rec))
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
        # Store the first number on the size list as the first of the incremental list
        self.incremental_sizes = []
        self.incremental_sizes.append(self.sizes[0])
        self.NG50 = 0
        self.LG50 = 0
        self.N50 = 0
        self.L50 = 0
        self.N90 = 0
        self.L90 = 0
        self.N95 = 0
        self.L95 = 0
        self.counter_over_1000 = 1 if (self.incremental_sizes[-1] > 1000) else 0
        # Now iterate over the sizes list, and incrementally add it to the incremental list
        for index, size in enumerate(islice(self.sizes, 1, None)):
            self.incremental_sizes.append(self.incremental_sizes[-1] + size)
            # While we iterate over the array, we can gather NG50, N50, LG50 and the L50 when we reach half of the total assembly
            if (expected_genome_size != 0):
                if (self.incremental_sizes[-1] > expected_genome_size / 2) and (self.NG50 == 0):
                    self.NG50 = size
                    # We need to add one to the index because we started to loop at 1, not 0
                    self.LG50 = index + 1
            else:
                self.NG50 = 0
                self.LG50 = 0
            if (self.incremental_sizes[-1] > self.total_length / 2) and (self.N50 == 0):
                self.N50 = size
                # We need to add one to the index because we started to loop at 1, not 0
                self.L50 = index + 1
            if (self.incremental_sizes[-1] > self.total_length * 0.90) and (self.N90 == 0):
                self.N90 = size
                # We need to add one to the index because we started to loop at 1, not 0
                self.L90 = index + 1
            if (self.incremental_sizes[-1] > self.total_length * 0.95) and (self.N95 == 0):
                self.N95 = size
                # We need to add one to the index because we started to loop at 1, not 0
                self.L95 = index + 1
            # Count the number of sequences larger than 1000
            if (size >= 1000):
                self.counter_over_1000 += 1

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
        return line

    def print_stats(self):
        # For prettyprinting the numbers
        import locale
        locale.setlocale(locale.LC_ALL, '')  # empty string for platform's default setting
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
        print "Name:%s" % self.name
        print "Count:%s" % format(self.count, "n")
        print "Sum:%s" % format(self.sum, "n")
        print "Max:%s" % format(self.max, "n")
        print "Min:%s" % format(self.min, "n")
        print "Average:%s" % format(self.average, "n")
        print "Median:%s" % format(self.median, "n")
        print "N50:%s" % format(self.N50, "n")
        print "L50:%s" % format(self.L50, "n")
        print "NG50:%s" % format(self.NG50, "n")
        print "LG50:%s" % format(self.LG50, "n")
        print "N90:%s" % format(self.N90, "n")
        print "L90:%s" % format(self.L90, "n")
        print "N95:%s" % format(self.N95, "n")
        print "L95:%s" % format(self.L95, "n")
        print "Count > 1000:%s" % format(self.counter_over_1000, "n")

# Check what kind of input we get
# Starts with ">" and is a single entry : 1 file
# Starts with ">" and there are more entries : multiple files
# Doesn't start with ">" : csv with 1 or more entries.
if (len(sys.argv) == 1):
    sys.exit()
if (len(sys.argv) == 2):
    # Single file, now check if the input is a FASTA formatted file
    with open(sys.argv[1], 'r') as f:
        first_line = f.readline()
    if (first_line[0] == '>'):  # We got a single FASTA file
        input_file = sys.argv[1]
        assemblies[input_file] = Assembly(input_file, input_file,min_length)
    else:  # We assume a valid CSV file, put the values in the assemblies dict
        file_name = sys.argv[1]
	fp = open(file_name, 'rb')
        csv_file = csv.DictReader((row for row in fp if not row.startswith('#')), delimiter=',', quotechar='"')
        for row in csv_file:
            input_name = row['Seq_Name']
            input_file = row['Seq_File']
            assemblies[input_name] = Assembly(input_name, input_file,min_length)
if (len(sys.argv) > 2):
    # Multiple files, for now assume they are all valid FASTA
    for input_file in sys.argv[1:]:
        assemblies[input_file] = Assembly(input_file, input_file,min_length)

# Now plot the A50 plots and print the stats
print csv2string(["Name","Count","Sum","Max","Min","Average","Median","N50","L50","NG50","LG50","N90","L90","N95","L95","Count > 1000"])
for name, assembly in iter(sorted(assemblies.items())):
    print(csv2string(assembly.return_stats()))
    line = pylab.plot(assembly.incremental_sizes, label=name,color=next(colors))

pylab.title("A50 plot of contigs >"+str(min_length)+"bp")
if (TITLE != '') : pylab.suptitle(TITLE, fontsize=20)
pylab.xlabel("Sequence count")
pylab.ylabel("Incremental size (bp)")
pylab.legend(loc='best')
if (max_contigs > 0):
    pylab.xlim([0, max_contigs])
if (outfile != ''):
    pylab.savefig(outfile,dpi = (DPI))
else:
    pylab.show()
