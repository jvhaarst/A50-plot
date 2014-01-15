#!/usr/bin/env python
# This scripts calculates some genome statistitics, and generates a so called A50 plot.
# Needed modules
from __future__ import division
import sys
from Bio import SeqIO
import pylab
from itertools import islice
# For prettyprinting the numbers
import locale
locale.setlocale(locale.LC_ALL, '') # empty string for platform's default setting
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

# Read in genome, and store lenghts in list
# TODO : make this a command line switch, and check for:
# Starts with ">" and is a single entry : 1 file
# Starts with ">" and there are more entries : multiple files
# Doesn't start with ">" : csv with 1 or more entries.

input_file = "/Users/jvhaarst/Dropbox/code/A50-plot/usda/usda/newbler-pbjelly.fasta"
input_name = "newbler+pbjelly"
sys.stderr.write("Reading in %s\n" % input_file)
# TODO : out from here till just before plotting into class
# For each contig/scaffold, calculate length, and sort large to small
sizes = [len(rec) for rec in SeqIO.parse(open(input_file), "fasta")]
sizes.sort(reverse=True)
# Calculate all the statistics we can already do
total_length = sum(sizes)
print "Count: %s"	%	format(len(sizes), "n")
print "Sum:%s"		%	format(total_length, "n")
print "Max:%s"		%	format(max(sizes), "n")
print "Min:%s"		%	format(min(sizes), "n")
print "Average:%s"	%	format(sum(sizes)/len(sizes), "n")
print "Median:%d" 	%	sizes[int(len(sizes)/2)]
# Store the first number on the size list as the first of the incremental list
incremental_sizes=[]
incremental_sizes.append(sizes[0])
N50=0
counter_over_1000=1 if (incremental_sizes[-1] > 1000) else 0
# Now iterate over the sizes list, and incrementally add it to the incremental list
for index,size in enumerate(islice(sizes, 1, None)):
	incremental_sizes.append(incremental_sizes[-1]+size)
	# While we iterate over the array, we can gather N50, L50 when we reach half of the total assembly
	if (incremental_sizes[-1] > total_length /2) and (N50==0):
		print "N50:%s" % format(size, "n")
		# We need to add one to the index because we started to loop at 1, not 0
		print "L50:%s" % format(index+1, "n")
		N50=size
	# Count the number of sequences larger than 1000
	if (size >=1000):
		counter_over_1000 += 1
print "Count > 1000:%s" % format(counter_over_1000, "n")

# Now plot the A50 plots.
pylab.plot(incremental_sizes, label=input_name)
pylab.title("A50 plot" )
pylab.xlabel("Sequence count")
pylab.ylabel("Incremental size (bp)")
pylab.legend(loc='best')
pylab.show()