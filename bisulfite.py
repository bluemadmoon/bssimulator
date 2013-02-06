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
from Bio.Alphabet import IUPAC



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


def bs_transform(seq, Methyl_positions, read_position, read_length):
# Bisulfite transformation function
# Transform all unmethylated Cs to Ts


    # Empty list which will contain the transformed sequence
    res_seq = list(seq.replace('C', 'T'))



    # For every C (i = position, w = methylation probability) ...
    for i, w in Methyl_positions:
        # If there is a methyl in the current read ...
        if i in xrange(read_position, read_position + read_length):
            if w_choice(w):
               res_seq[int(i - read_position)] = 'C'

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
    
    position = np.vstack((np.array(position), np.random.rand(len(position))))
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

    # Argument to set the length of the reads to be sequenced
    parser.add_option("-R", 
        "--read_length", 
        action="store", 
        help="Size of the reads (must be >10 and <100, set to 50 by default)", 
        default=50)


    parser.add_option("-S",
        "--seq", 
        action="store", 
        help="Number of sequencing runs (higher number gives better coverage, but much \
            more computation time (not recommended to go above 200, default = 100) ", 
        default=100)

    (options, args) = parser.parse_args() 

    if not options.input:
        sys.stdout.write("Sorry: you must specify an input file\n")
        sys.stdout.write("More help avalaible with -h or --help option\n")
        sys.exit(0)

    if options.read_length < 10 or options.read_length > 100:
        sys.stdout.write("Read length incorrect.\n Must be in [10, 100] range\nQuitting")
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

        print " - Computing random methylation"
        island_position = IslandPosition(length, [100, 800], [1000, 3000])
        Methyl_pos = MethylPosition(island_position)

        # A not so dirty to write data to a file !
        for i in xrange(len(Methyl_pos)):
            profile_output_file.write("%i\t%.3f\n" % (Methyl_pos[i, 0], Methyl_pos[i, 1]))
        profile_output_file.close()

        print " - Computing the bisulfite process and sequencing"

        readLength = options.read_length
        nbSeq = int(options.seq)
        print " - Simulating reads sequencing"
        print "   Read length : %d" % int(readLength)
        print "   Number of sequencings : %d" % int(nbSeq)
        k = 0
        for i in xrange(nbSeq): # Number of sequencings.
            read_positions = read_generation(length, int(0.05*length), readLength)
            for i in read_positions:

                cur_read = cur_record.seq[i:i + readLength]
                bs_read = bs_transform(str(cur_read), Methyl_pos, i, readLength)

                output_file.write(">%s|r%i\n"%(cur_record.name, int(k)))
                k = k + 1
                # we take only the 50 elements from each reads
                output_file.write(str(bs_read)+'\n')
    
    input_file.close()
    output_file.close()
    print "Process began at %s : " % begin_time
    print "Process terminated at %s : " % strftime("%Y-%m-%d %H:%M:%S")
    
if __name__ == '__main__':
    main()