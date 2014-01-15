#!/usr/bin/env python -t
# This scripts calculates some genome statistitics, and generates a so called A50 plot.
# Needed modules
from __future__ import division
import csv
import sys
import pylab
from Bio import SeqIO
from itertools import islice
# For prettyprinting the numbers
import locale
locale.setlocale(locale.LC_ALL, '') # empty string for platform's default setting
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
# Global variables
assemblies = {}
outfile='assemblies.png' # Leave empty if you want to use the interactive output instead of a file.

# Class definitions
class Assembly:
	"""Class to hold the assembly information, and perform the calculations"""
	def __init__(self, name, filename):
		self.name = name
		self.filename = filename
		sys.stderr.write("Reading in %s\n" % self.filename)
		# For each contig/scaffold, calculate length, and sort large to small
		self.sizes = [len(rec) for rec in SeqIO.parse(open(self.filename ), "fasta")]
		self.sizes.sort(reverse=True)
		# Calculate all the statistics we can already do
		self.total_length	= sum(self.sizes)
		self.count 			= len(self.sizes)
		self.sum			= self.total_length
		self.max 			= max(self.sizes)
		self.min 			= min(self.sizes)
		self.average 		= sum(self.sizes)/len(self.sizes)
		self.median 		= self.sizes[int(len(self.sizes)/2)]
		# Store the first number on the size list as the first of the incremental list
		self.incremental_sizes=[]
		self.incremental_sizes.append(self.sizes[0])
		self.N50=0
		self.L50=0
		self.counter_over_1000=1 if (self.incremental_sizes[-1] > 1000) else 0
		# Now iterate over the sizes list, and incrementally add it to the incremental list
		for index,size in enumerate(islice(self.sizes, 1, None)):
			self.incremental_sizes.append(self.incremental_sizes[-1]+size)
			# While we iterate over the array, we can gather N50, L50 when we reach half of the total assembly
			if (self.incremental_sizes[-1] > self.total_length /2) and (self.N50==0):
				self.N50=size
				# We need to add one to the index because we started to loop at 1, not 0
				self.L50=index+1
			# Count the number of sequences larger than 1000
			if (size >=1000):
				self.counter_over_1000 += 1

# Check for what kind of input we get
# Starts with ">" and is a single entry : 1 file
# Starts with ">" and there are more entries : multiple files
# Doesn't start with ">" : csv with 1 or more entries.
if (len(sys.argv) == 1):
	sys.exit()
if (len(sys.argv) == 2):
	# Single file, now check if the input is a FASTA formatted file
	with open(sys.argv[1], 'r') as f:
		first_line = f.readline()
	if (first_line[0] == '>'): # We got a single FASTA file
		input_file = sys.argv[1]
		assemblies[input_file] = Assembly(input_file,input_file)
	else: # We assume a valid CSV file, put the values in the assemblies dict
		file_name = sys.argv[1]
		csv_file = csv.DictReader(open(file_name, 'rb'), delimiter=',', quotechar='"')
		for row in csv_file:
			input_name = row['Seq_Name']
			input_file = row['Seq_File']
			assemblies[input_name] = Assembly(input_name,input_file)
if (len(sys.argv) > 2):
	# Multiple files, for now assume they are all valid FASTA
	for input_file in sys.argv[1:]:
		assemblies[input_file] = Assembly(input_file,input_file)

# TODO : print the statistics as below for each Assembly object.
# print "Count: %s"	%	format(len(sizes), "n")
# print "Sum:%s"		%	format(total_length, "n")
# print "Max:%s"		%	format(max(sizes), "n")
# print "Min:%s"		%	format(min(sizes), "n")
# print "Average:%s"	%	format(sum(sizes)/len(sizes), "n")
# print "Median:%d" 	%	sizes[int(len(sizes)/2)]
# print "N50:%s" % format(size, "n")
# print "L50:%s" % format(index+1, "n")
# print "Count > 1000:%s" % format(counter_over_1000, "n")

# Now plot the A50 plots.
for name, assembly in assemblies.iteritems():
	pylab.plot(assembly.incremental_sizes, label=name)
pylab.title("A50 plot" )
pylab.xlabel("Sequence count")
pylab.ylabel("Incremental size (bp)")
pylab.legend(loc='best')
if (outfile != ''):
	pylab.savefig(outfile)
else:
	pylab.show()