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

Use -h (or --help) argument for further help on usage and possible arguments .
"""

import sys, os, random, re
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


def w_choice(prob):
# Function for boolean random weight picking (True : methylated C)


    lst = [(True, prob), (False, (1 - prob))]
    n = random.uniform(0, 1)
    for state, weight in lst:
        if n < weight:
            break
        n = n - weight
    return state



def bs_transform(seq, Methyl_positions, read_length):
# Bisulfite transformation function
# Transform all unmethylated Cs to Ts


    # Empty list which will contain the transformed sequence
    res_seq = list(seq.replace('C', 'T'))

    # For every C (i = position, w = methylation probability) ...
    for i, w in Methyl_positions:
        # If there is a methyl in the current read ...
        if w_choice(w):
            res_seq[int(i)] = 'C'

    # Return the transformed sequence as str
    return ''.join(res_seq)


def IslandPosition(length, taille_moy, distance_moy):

    tab = np.zeros((2, 2))  # initialisation [position_centre, longueur_total]

    # tant que la derniere position + la taille est inferieure a la taille totale
    while int(tab[-1][0] + tab[-1][1] + distance_moy[1]) < length:

        # on tire une distance entre les centres
        distance = np.random.randint(distance_moy[0], distance_moy[1])

        # on calcule la nouvelle distance
        position = tab[-1][0] + distance
         # on tire une taille
        taille = np.random.randint(taille_moy[0], taille_moy[1])

        # on checke les conditions
        # si la nouvelle borne inferieure est superieur a la derniere borne superieur
        if int(position - (taille / 2)) > int(tab[-1][0] + (tab[-1][1] / 2) + 1):
            tab = np.vstack((tab, [position, taille]))  # on ajoute la nouvelle position

    tab = np.delete(tab, [0, 1], 0)
    return tab

def MethLevelSample(size):
    # Choose a different level of methylation weighted by a given probability
    # Using the same weights as in Martin C. et al. article (A mostly
    # traditional approach improves alignemnt of bisulfite-converted DNA)
    # We assign 5 distinct levels of methylation, and we assume that we have 
    # higher methylation probability in CpG islands, thus allowing for some 
    # Cs to be methylated with lower probability

    sample = []
    for i in range(size):
        weights = [.05, .05, .05, .05, .8]
        levels = [0, .1, .2, .5, 1]
        sample.append(np.random.choice(levels, p=weights))
    return sample

def MethylPosition(tab):
# Compute random methylated cytosine position


    # retourne la position des methyl dans le genome en suivant la position des island
    position = []
    for i,j in tab:
        CGindex = np.array([m.start() for m in re.finditer \
            ('CG', str(cur_record.seq[int(i - j / 2): int(i + j / 2)]))])
        if len(CGindex) <> 0:  # conditions to check if there are cpg in the island
            a, b = np.histogram(np.random.random(len(CGindex)), len(CGindex))
            position.extend(CGindex[np.where(a>0)]+(i-(j/2)))
    
    position = np.vstack((np.array(position), MethLevelSample(len(position))))
    return position.T

def read_generation(gen_len, rest_nb, read_length):
# 100 to 200bp reads generation by single strand DNA cutting by 
# simulating restriction enzymes

    # Randomly choose rest_nb restriction sites along the sequence
    # We sort them and make them unique
    cuts = np.unique(np.sort(np.random.randint(gen_len, size = rest_nb)))


    # We only keep reads of size ranging from 100 to 200 bp
    # Also, trim the reads to make the comply to the selected size by the user
    # This is not very useful in the case of a simulation, but in real data,
    # reads are often trimmed car the end of the reads tends to have lower
    # sequençing quality, besides the need to have equal size reads.
    reads = []
    for i in xrange(1, len(cuts)):
          if (cuts[i] - cuts[i - 1]) > 100 and (cuts[i] - cuts[i - 1]) < 200:
            reads.append(cuts[i])

    return reads
    
    
def main():
    # -----------------------------------
    # ARGUMENT MANAGER
    # -----------------------------------
    
    # Usual verifications and warnings (errors)
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

    # Argument to set the length of the reads to be sequenced
    parser.add_option("-R", 
        "--read_length", 
        action="store", 
        type=int,
        help="Size of the reads (must be >10 and <100, set to 50 by default)", 
        default=50)

    # Argument to set the number of sequencing runs
    parser.add_option("-S",
        "--seq", 
        action="store", 
        type=int,
        help="Number of sequencing runs (higher number gives better coverage, but much \
            more computation time (not recommended to go above 200, default = 100) ", 
        default=100)

    (options, args) = parser.parse_args() 

    if not options.input:
        sys.stdout.write("Sorry: you must specify an input file\n")
        sys.stdout.write("More help avalaible with -h or --help option\n")
        sys.exit(0)

    if options.read_length < 10 or options.read_length > 100:
        sys.stdout.write("Read length incorrect.\n Must be in [10, 100] range\nQuitting\n")
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
    output_filename = os.getcwd() \
    + '/' \
    + os.path.splitext(base)[0] \
    + '_bs.fasta' 

    profile_output_filename = os.getcwd() \
    + '/' \
    + os.path.splitext(base)[0] \
    + '_profile.dat'    
    
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

        print " - Sampling random methylation"
        island_position = IslandPosition(length, [100, 800], [1000, 3000])
        methyl_pos = MethylPosition(island_position)

        # Saving the methylation profile.
        # This is the actual methylation state of the DNA, the actual information
        # we will try to infer. We're keeping it for validation purpose 
        # A not so dirty to write data to a file !
        for i in xrange(len(methyl_pos)):
            profile_output_file.write("%i\t%.3f\n" % (methyl_pos[i, 0], methyl_pos[i, 1]))
        profile_output_file.close()

        print " - Computing the bisulfite process and sequencing"

        read_length = options.read_length
        nbSeq = options.seq

        print "   Read length : %d" % int(read_length)

        k = 0
        for i in xrange(nbSeq): # Number of sequencings.

            sys.stdout.write("\r   Sequencing progress : %i/%i" %(i, nbSeq)) 
            sys.stdout.flush()

            bs_record = Seq(bs_transform(str(cur_record.seq),
                methyl_pos,
                read_length))

            read_positions = read_generation(length, int(0.01*length), read_length)
            for j in read_positions:

                cur_read = bs_record[j:j + read_length]

                output_file.write(">%s|r%i\n"%(cur_record.name, k))
            
                output_file.write(str(cur_read)+'\n')
                k += 1

    input_file.close()
    output_file.close()
    print "\nProcess began at %s : " % begin_time
    print "Process terminated at %s : " % strftime("%Y-%m-%d %H:%M:%S")
    
if __name__ == '__main__':
    main()