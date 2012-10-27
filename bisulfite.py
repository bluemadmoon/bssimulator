#!/usr/bin/python
# encoding: utf-8
"""
bisulfite.py

Simulate a sodium bisulphite transformation of a single DNA strand (without PCR).
Takes a fasta file as input, generates random methylated cytosines
and outputs a files containing reads, ready to be aligned with
a modified index (using bisulfite-build.py)

Usage :
$python bisulfite-build.py -i [input file]

Use -h argument for further help on usage.

Copyright (c) 2012 Madjid BESSOUL and Guillaume VIEJO. All rights reserved.
"""

import sys
import os
from optparse import OptionParser
from time import strftime
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

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Tryin to import Matplotlib
try:
    import matplotlib.pyplot as plt
    mpl = True
except ImportError:
    mpl = False
    print "Matplotlib package is absent, you can't visualise histograms"

import random
import re


# Function for boolean random weight picking (True : methylated C)
def w_choice(prob):

    lst = [(True, prob), (False, (1 - prob))]
    n = random.uniform(0, 1)
    for state, weight in lst:
        if n < weight:
            break
        n = n - weight
    return state


def query_yes_no(question, default="yes"):

    valid = {"yes": True,   "y": True,  "ye": True,
             "no": False,     "n": False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


def bs_transform(seq, pos):
    # Bisulfite transformation function

    # Empty list which will contain the transformed sequence
    res_seq = list(seq.replace('C', 'T'))

    # For every C (i = position, w = methylation probability) ...
    for i, w in pos:
        if w_choice(w):
            res_seq[int(i)] = 'C'

    # Return the transformed sequence as str
    return ''.join(res_seq)


def IslandPosition(length, taille_moy, distance_moy):
    # entree : taille total , taille moyenne = [min, max], distance moyenne = [min, max],
    # then tirage

    tab = np.zeros((2, 2))  # initialisation [position_centre, longueur_total]

    # tant que la derniere position + la taille est inferieure a la taille totale
    while int(tab[-1][0] + tab[-1][1] + distance_moy[1]) < length:

        # on tire une distance entre les centres
        distance = np.random.randint(distance_moy[0], distance_moy[1])

        # on calcule la nouvelle distance
        position = tab[-1][0]+distance
        taille = np.random.randint(taille_moy[0], taille_moy[1]) # on tire une taille
        
        # on checke les conditions
        # si la nouvelle borne inferieure est superieur a la derniere borne superieur   
        if int(position-(taille/2)) > int(tab[-1][0]+(tab[-1][1]/2)+1):
            tab = np.vstack((tab, [position, taille])) # on ajoute la nouvelle position
    
    tab = np.delete(tab, [0,1],0)
    return tab

# Compute random methylated cytosine position
def MethylPosition(tab):
    # retourne la position des methyl dans le genome en suivant la position des island
    position = []
    for i,j in tab:
        CGindex = np.array([m.start() for m in re.finditer('CG', str(cur_record.seq[int(i-j/2):int(i+j/2)]))])
        if len(CGindex) <> 0:  # conditions to check if there are cpg in the island
            a, b = np.histogram(np.random.random(len(CGindex)), len(CGindex))
            position.extend(CGindex[np.where(a>0)]+(i-(j/2)))
    
    position = np.vstack((np.array(position), np.random.rand(len(position))))
    return position.T

# Fuction that simulates a 50bp long read generation
def read_generation(L, N):
    # input: length of the genome, number of restriction site
    # output: [position min, position max] of read to extract in the genome
    # Not all the genome is extracted since only reads ranging from 50 to 100 bp are selected.
    site = np.sort(np.random.randint(L, size = N))
    #taille = np.array([site[i]-site[i-1] for i in xrange(1, len(site))])
    site2 = []
    # selection of only 50-100 bp
    for i in xrange(1, len(site)):
          if (site[i]-site[i-1])>50 and (site[i]-site[i-1])<100:
            site2.append([site[i-1], site[i]])
              
    return np.array(site2)
    
    
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
    parser.add_option("-i", 
        "--input", 
        action="store", 
        help="The genome input file, must be a fasta or fastq file", 
        default=False)


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
    except IOError :
        print "No such file found. \nClosing."
        sys.exit(0)

    input_file = open(options.input, 'r') 
    
    # Retrive input filename base
    base = os.path.basename(options.input)
    
    # Create an input_file_index.fasta name
    output_filename = os.getcwd() + '/' + os.path.splitext(base)[0] + '_bs.fasta'   
    profile_output_filename = os.getcwd() + '/' + os.path.splitext(base)[0] + '_profile.dat'    
    
    # File to store the reads
    output_file = open(output_filename,'a')
    
    # File to store the profile
    profile_output_file = open(profile_output_filename, 'a')

    global cur_record
    for cur_record in SeqIO.parse(input_file, "fasta"):

        length = len(cur_record.seq)

        # Store the begin time
        begin_time = strftime("%Y-%m-%d %H:%M:%S")

        print "Processing gene : %s" % cur_record.name

        print " - Computing random methylation"
        island_position = IslandPosition(length, [100, 800], [1000, 3000])
        pos = MethylPosition(island_position)

        # A not so dirty to write data to a file !
        for i in xrange(len(pos)):
            profile_output_file.write("%i\t%.3f\n" % (pos[i, 0], pos[i, 1]))
        profile_output_file.close()

        print " - Computing the bisulfite process"
        cur_record.seq = Seq(bs_transform(str(cur_record.seq), pos), IUPAC.unambiguous_dna)
        print " - Simulating reads sequencing"
        k = 0
        for i in xrange(100): # Number of sequencings.
            site = read_generation(length, int(0.05*length))
            for i in site[:,0]:
                  output_file.write(">%s|r%i\n"%(cur_record.name, int(k)))
                  k = k+1
                  # we take only the 50 elements from each reads
                  output_file.write(str(cur_record.seq[i:i+50])+'\n')
    
    input_file.close()
    output_file.close()
    print "Process began at %s : " % begin_time
    print "Process terminated at %s : " % strftime("%Y-%m-%d %H:%M:%S")
    
if __name__ == '__main__':
    main()