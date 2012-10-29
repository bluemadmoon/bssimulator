#!/usr/bin/python
# encoding: utf-8
"""
bisulfite-build.py

Build new index for BS data alignment using Bowtie
Takes a fasta file as input.
Computes a C to T transformation and returns a fasta file.

Usage :
$python bisulfite-build.py -i [input file]

Use -h argument for further help on usage.



Copyright (c) 2012 Madjid BESSOUL and Guillaume VIEJO. All rights reserved.
"""



import sys
import os
from optparse import OptionParser

# Trying to import Numpy package # REQUIRED
try:
	import numpy as np
except ImportError:
	print "Numpy package is required to launch the program. \nExitting"
	sys.exit(0)
	
# Trying to import Biopython package # REQUIRED	
try:		
	from Bio import SeqIO
except ImportError:
	print "Biopython package is required to launch the program. \nExitting"
	sys.exit(0)
# Tryin to import Matplotlib
try:		
	import matplotlib.pyplot as plt
	mpl = True
except ImportError:
	mpl = False
	print "Matplotlib package is absent, you can't visualise histograms"

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import random
import re


def main():
	
	# -----------------------------------
	# ARGUMENT MANAGER
	# -----------------------------------
	
	# Usual verifications and warnings
	if not sys.argv[1:]:
		sys.stdout.write("Sorry: you must specify at least 1 argument")
		sys.stdout.write("More help avalaible with -h or --help option")
		sys.exit(0)

	parser = OptionParser()
	
	# Argument to fetch the input filename.
	parser.add_option("-i", "--input", action="store", help="The genome input file, must be a fasta or fastq file", default=False)
			
	
	parser.add_option("-s", "--stats", action="store_true", help="Show histograms and statistics", default=False )

	(options, args) = parser.parse_args() 


	if not options.input:
		sys.stdout.write("Sorry: you must specify an input file\n")
		sys.stdout.write("More help avalaible with -h or --help option\n")
		sys.exit(0)
		
	# -----------------------------------
	# END ARGUMENT MANAGER
	# -----------------------------------

	# File existence test
	try:
		open(options.input, 'r')
	except:
		print "No such file found. \nClosing."
		sys.exit(0)
		
	input_file = open(options.input, 'r')
	
	# Retrive input filename base
	base = os.path.basename(options.input)
	# Create an input_file_index.fasta name
	output_filename = os.getcwd() + '/' + os.path.splitext(base)[0] + '_index.fasta'	
	output_file = open(output_filename,'w')
	
	
	for cur_record in SeqIO.parse(input_file, "fasta"):
				
		### Show gene statistics
		if options.stats and mpl:
			gene_name = cur_record.name 
			C_count = cur_record.seq.count('C') 
			G_count = cur_record.seq.count('G') 
			seq_length = len(cur_record.seq) 
			cg_percentage = float(C_count + G_count) / seq_length 
		
			# How much CpG through the whole genome
			CGcontent = cur_record.seq.count('CG')
		
			# Get the position of every single CpG
			GCindex = [m.start() for m in re.finditer('CG', str(cur_record.seq))]
	
			# Parameters
			length = 50000
			threshold = 400 # Show only steps with more than 10 CpGs
			step = 1000	# length of the histogram step

			# Histogram
			fig = plt.figure()
			plt.hist(GCindex, bins=step, normed=False, facecolor="green")
		
			ax = [0, seq_length, threshold, 1000]
			plt.axis(ax)
		
			plt.grid(True)
			txt = 'Gene name : %s\nGCpercentage = %f\nNumber of CpG locations = %i ' % (gene_name, 	cg_percentage, len(GCindex))
			fig.text(.15, .75, txt)
		
			plt.title(r'CpG concentration, step = %i bp' % step)
			plt.show()
			
		# Create a new index with only A, T and G and store 
				
		cur_record.seq = Seq(str(cur_record.seq).replace('C', 'T'), IUPAC.unambiguous_dna)
		SeqIO.write(cur_record, output_file, "fasta")

	output_file.close()
	input_file.close()
if __name__ == '__main__':
    main()