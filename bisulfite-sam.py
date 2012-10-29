#!/usr/bin/python
# encoding: utf-8
"""
bisulfite-sam.py
Analysing the sam alignment data
For each C position, find the number of misaligned nucleotides and the number of 
overlapping reads at this position.

The BAM/SAM file must be indexed.

/!\ REQUIRED PACKAGES : 
		- Pysam : SAM/BAM file parser
		- Numpy : Numerical analysis package
Usage : 

python bisulfite-sam.py -i [bam_file] -r [fasta_file]

Use the option -h for help.
Copyright (c) 2012 Madjid BESSOUL and Guillaume VIEJO. All rights reserved.
"""

import sys
import os
import time
import datetime
from optparse import OptionParser
	
# Trying to import Biopython package # REQUIRED	
try:		
	from Bio import SeqIO
except ImportError:
	print "Biopython package is required to launch the program. \nExitting"
	sys.exit(0)
from Bio.Seq import Seq

try:
	import pysam as ps
except ImportError:
	print "Pysam package is required to launch the program. \nExitting"
	sys.exit(0)
	
try:
	import numpy as np
except ImportError:
	print "Numpy package is required to launch the program. \nExitting"
	sys.exit(0)
	
# Splitting function so we can parse the MD tag (faster than regex parsing)
# We need to retrieve the numbers seperated by ATGCN
def split(tag, seps):
    default_sep = seps[0]
    for sep in seps[1:]: # we skip seps[0] because that's the default seperator
        tag = tag.replace(sep, default_sep)
    return [i.strip() for i in tag.split(default_sep)]

def main():
	#-----------------------------------
	# ARGUMENT MANAGER
	# -----------------------------------
	
	# Usual verifications and warnings
	if not sys.argv[1:]:
		sys.stdout.write("Sorry: you must specify at least 1 argument")
		sys.stdout.write("More help avalaible with -h or --help option")
		sys.exit(0)
	
	parser = OptionParser()
	
	# Argument to fetch the input filename.
	parser.add_option("-i",
	                 "--input", 
	                 action="store",
	                 help="The input file must be a BAM/SAM file", 	 
	                 default=False)
	
	parser.add_option("-r",
	 				"--ref", 
					action="store",
					help="The reference genome file must be a FASTA/FASTQ file",
					default=False)
	
	(options, args) = parser.parse_args() 
	
	if not options.input:
		sys.stdout.write("Sorry: you must specify an input file\n")
		sys.stdout.write("More help avalaible with -h or --help option\n")
		sys.exit(0)
		
	# -----------------------------------
	# END ARGUMENT MANAGER
	# -----------------------------------
	
	print """
	bisulfite-sam.py
	
	"""
	
	# -----------------------------------
	# END WELCOME MESSAGE
	# -----------------------------------
	overall_start = time.clock()
	samfile = ps.Samfile(options.input, 'rb') 
	
	base = os.path.basename(options.input)
	output_filename = os.getcwd() + '/' + os.path.splitext(base)[0] + '_meth_pos.dat'
	outfile = open(output_filename,'w')
	
	
	
		# Opening the fasta file
	ref_file = open(options.ref, 'r')
	
	# We find the cytosines in the genomes, and output the mismatches on these
	# positions only.
	
	print " - Building cytosine positions index"
	
	cyt = []
	for cur_record in SeqIO.parse(ref_file, "fasta"):
		for i in xrange(len(cur_record.seq)):
			if cur_record.seq[i] == 'C':
				cyt.append(i)
	
	# Var init
	mismatch = {} 
	coverage = {}

	print ' - Counting mismatches and overlapping reads over BAM file'
	start = time.clock()
	# Main loop, reading the sam file line by line (ie. read by read)
	for read in samfile.fetch():
		### Coverage counting
		if read.pos not in coverage:
			coverage.update({read.pos:1})
		else:
			coverage[read.pos] += 1
		
		if read.tags[1][1] <> '50':
			single_mismatch_pos = 0
			for i in split(read.tags[1][1], ['T', 'G', 'A', 'C', 'N'])[0:-1]:
			
				single_mismatch_pos += int(i) 
				pos = read.pos + single_mismatch_pos 
				
				# New position of mismatch
				if pos in mismatch: 
					# update number of mismatch
					mismatch[pos] += 1
				else :
					# else create one 
					mismatch[pos] = 1
					
	
	end = time.clock()
	print '     done in ', end - start, 'seconds.'				
	print ' - Building data and writing to file'
	start = time.clock()
	
	
	# Add two columns to the cyt array, to contain converage and mismatches
	# Cyt will contain the final data
	cyt = np.column_stack([cyt, [0] * len(cyt)])
	cyt = np.column_stack([cyt, [0] * len(cyt)])

	# Sorting the coverage keys, hoping it can speed up things a bit ...
	keylist = sorted(coverage.keys())
	
	# Write header to file
	
	outfile.write(u"# bisulfite-sam.py output for %s\n" % base)
	outfile.write("# %s\n" % datetime.datetime.now().ctime())
	outfile.write(u"""
# column 1 : position of every cytosine over the genome
# column 2 : overlapping at the cyt position (number of overlapping reads)
# column 3 : number of mismatches by bowtie alignment over C-less DNA index

""")

	
	for i in xrange(len(cyt)):
		
		# Count mismatches at the C positons ...
		if cyt[i, 0] in mismatch:
			cyt[i, 2] = mismatch[cyt[i, 0]]
			
		# /!\ todo ...	
		# And the coverage ... still getting weird results. 
		for key in keylist:
			if key < cyt[i, 0] and key >= cyt[i,0] - 50:
				cyt[i, 1] += coverage[key]
				
			# Break if out of bounds
			if key > cyt[i, 0]:
				break
				
		# Writing it down line by line baby !

			
		outfile.write("%i\t\t%i\t\t%i\n" % (cyt[i, 0], cyt[i, 1], cyt[i, 2]))			
		
   
	end = time.clock()
	print '     done in', end - start, 'seconds.'
	overall_end = time.clock()
	print 'Overall process run time : ', overall_end - overall_start
	print 'Quitting.'
	
	
	samfile.close()
	outfile.close()
	ref_file.close()
if __name__ == '__main__':
	main()