#!/usr/bin/env python3
# Name: Sherry Lin
#
import sys
import time
from random import seed
from random import randint


"""
Homework 2: Finding CRISPR arrays

Input/Output: STDIN / STDOUT


Examples:


"""

# TODO Implement RanonizedMotif Search
# TODO Implement Gibbs sampling

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.

        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Parse arguments for search for the missing',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        # Main Assignment
        self.parser.add_argument('-i', '--iterations', type=int, help='The number of iterations for running the search')
        self.parser.add_argument('-k', '--motifLength', type=int, help='The length of the target motif')
        self.parser.add_argument('-p', '--pseudocount', type=float, help='Pseudocount when computing Profile')
        # Extra credit
        self.parser.add_argument('-r', '--shuffle', action='store_true', help='Shuffling the input sequence ')
        self.parser.add_argument('-g', '--gibbs', action='store_true', help='using gibbs sampling ')
        self.parser.add_argument('-m', '--matrix', action='store_true', help='print the name and the motif ')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class FastAreader():
    """
    Helper function that returns objects of header and the DNA sequences separately

    """

    def __init__(self, fname=''):
        self.fname = fname

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class SearchMotif():
    # TODO Randomly select one 5-mer motif from each sequence and generate a Motif Matrix
    # TODO Count nucleotides and generat a Profile
    # TODO Adding the pseudo count
    # TODO Compute Profile probability
    # TODO Compute profile entropy
    # TODO Generate new Motif Matrix by the highest probability motif in the sequence
    # TODO Count nucleotides and generat a Profile
    # TODO Adding the pseudo count
    # TODO Compute Profile probability
    # TODO Compute profile entropy
    # TODO Compare the two entropy, takes the sequence with lowest entropy, stop when entropy goes up
    def __init__(self, sequences, iterations, motifLength, pseudocount, shuffle, gibbs, matrix):
        self.sequences = sequences
        self.iterations = iterations
        self.motifLength = motifLength
        self.pseudocount = pseudocount
        self.shuffle = shuffle
        self.gibbs = gibbs
        self.matrix = matrix

    def selectRandomMotif(self):
        """
        Helper function to generate random Motifs with target motif length
        Returns:
            A motif matrix with motifs randomly selected from the DNA sequences

        """
        seed(914)
        motifs = []
        for sequence in self.sequences:
            # https://docs.python.org/3/library/random.html
            i = randint(0, len(sequence)-self.motifLength)
            # print(i)
            motifs.append(sequence[i:i+self.motifLength])
        return motifs

    def generateProfile(self, motifs):
        # https://stackoverflow.com/questions/4937491/matrix-transpose-in-python
        motifsTranspose = [*zip(*motifs)]
        AList, CList, GList, TList = [], [], [], []
        for col in motifsTranspose:
            AList.append(col.count('A') + self.pseudocount)
            CList.append(col.count('C') + self.pseudocount)
            GList.append(col.count('G') + self.pseudocount)
            TList.append(col.count('T') + self.pseudocount)
        print(AList, CList, GList, TList)

    def randomizedMotifSearch(self):
        pass

class Usage(Exception):
    """
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    """

    def __init__(self, msg):
        self.msg = msg


def main(myCommandLine=None):
    """
    Implement finding the missing sequence

    """

    try:
        myCommandLine = CommandLine()  # read options from the command line
        print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do
    except Usage as err:
        print(err.msg)

    # Get the commandline arguments
    iterations = myCommandLine.args.iterations
    motifLength = myCommandLine.args.motifLength
    pseudocount = myCommandLine.args.pseudocount

    shuffle = myCommandLine.args.shuffle
    gibbs = myCommandLine.args.gibbs
    matrix = myCommandLine.args.matrix

    print('extracted args', iterations, motifLength, pseudocount, shuffle, gibbs, matrix)

    fastaFile = FastAreader().readFasta()
    # store all sequence in a list
    sequences = []
    for header, sequence in fastaFile:
        #print('header is', header)
        #print('seq is', sequence)
        # print(len(sequence))
        sequences.append(sequence)
    print('DNA seqs are',sequences, 'length is', len(sequences))


if __name__ == "__main__":
    #start = time.time()
    main()
    #print('time consumed is', time.time() - start)
