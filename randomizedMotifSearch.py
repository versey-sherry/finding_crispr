#!/usr/bin/env python3
# Name: Sherry Lin
#
import sys
import time
from random import seed
from random import randint
from random import shuffle
from collections import defaultdict
from math import log2

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
        self.parser.add_argument('-r', '--scramble', action='store_true', help='Shuffling the input sequence ')
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
    def __init__(self, headers, sequences, iterations, motifLength, pseudocount, scramble, gibbs, matrix):
        self.headers = headers
        self.sequences = sequences
        self.iterations = iterations
        self.motifLength = motifLength
        self.pseudocount = pseudocount
        self.scramble = scramble
        self.gibbs = gibbs
        self.matrix = matrix

    def selectRandomMotif(self):
        """
        Helper function to generate random Motifs with target motif length
        Returns:
            A motif matrix with motifs randomly selected from the DNA sequences

        """
        # seed(914)
        motifs = []
        for sequence in self.sequences:
            # https://docs.python.org/3/library/random.html
            i = randint(0, len(sequence) - self.motifLength)
            # print(i)
            motifs.append(sequence[i:i + self.motifLength])
        return motifs

    def generateProfile(self, motifs):
        """
        Input motif matrix and generate a profile from the motif matrix.
        Args:
            motifs: The motif for generating the Profile

        Returns:
            A dictionary with the column number as key and nucleotides count dictionary as dictionary values

        """
        # the sum of nucleotide count per column including the added pseudocounts
        colSum = len(motifs) + 4 * self.pseudocount
        # print(colSum)
        # https://stackoverflow.com/questions/4937491/matrix-transpose-in-python
        motifsTranspose = [*zip(*motifs)]
        profileDict = defaultdict(dict)
        for colNum, col in enumerate(motifsTranspose):
            profileDict[colNum]['A'] = (col.count('A') + self.pseudocount) / colSum
            profileDict[colNum]['C'] = (col.count('C') + self.pseudocount) / colSum
            profileDict[colNum]['G'] = (col.count('G') + self.pseudocount) / colSum
            profileDict[colNum]['T'] = (col.count('T') + self.pseudocount) / colSum
        return profileDict

    def computeEntropy(self, profile):
        """
        Compute the entropy of the profile
        Args:
            profile: the profile of probability from the motif matirx

        Returns:
            the entropy of the profile

        """
        entropy = -sum([sum([probability * log2(probability) for probability in colValue.values()])
                        for colValue in profile.values()])
        return entropy

    def findNewMotifs(self, profile):
        """
        Given a profile, look at all motifs in each sequence and generate a new motif matrix
        with the most probable motif given profile
        Args:
            profile: the condition profile

        Returns:
            a new motif matrix with the most probable motif in each DNA sequence given profile

        """
        # print(profile)
        motifs = []
        for sequence in self.sequences:
            # initialize the best motif for the sequence
            pickMotif = "Z"
            maxProb = 0.
            for n in range(len(sequence) - self.motifLength + 1):
                tempSeq = sequence[n:n + self.motifLength]
                seqProb = 1.

                # computer probabilty for the particular motif
                for index, nucleotide in enumerate(tempSeq):
                    # print('index, nucleotide', index, nucleotide)
                    # print('prob', profile[index][nucleotide])
                    seqProb = seqProb * profile[index][nucleotide]

                # pick the motif with max probability and lower alpha order
                if seqProb >= maxProb:
                    maxProb = seqProb
                    pickMotif = tempSeq

            motifs.append(pickMotif)
        return motifs

    def randomizedMotifSearch(self):
        """
        Put everything together and implement randomizedMotifSearch
        Returns:
            bestMotifs: The motif matrix that's associated with the best score
            bestScore: The best entropy score associated with the bestMotifs

        """
        motifs = self.selectRandomMotif()
        # Initialize the best set of Motifs
        bestMotifs = motifs
        bestProfile = self.generateProfile(motifs)
        bestScore = self.computeEntropy(bestProfile)

        while True:
            tempMotifs = self.findNewMotifs(bestProfile)
            tempProfile = self.generateProfile(tempMotifs)
            tempScore = self.computeEntropy(tempProfile)
            # print('Best score and temp score', bestScore, tempScore)
            # Compare the scores
            if tempScore < bestScore:
                bestMotifs = tempMotifs
                bestProfile = tempProfile
                bestScore = tempScore
            else:
                return bestMotifs, bestScore

    def gibbsSampler(self):
        pass

    def iterateSearch(self):
        """
        Iterate the search a lot of times and find the best score
        Returns:
            the motif matrix with the best associated with the best score, the concensus motif and the best score
        """
        # Initialized bestScore with a big number
        allBestScore = 1000000000
        allBestMotifs = []

        for i in range(self.iterations):
            bestMotifs, bestScore = self.randomizedMotifSearch()
            if bestScore < allBestScore:
                allBestScore = bestScore
                allBestMotifs = bestMotifs
        consensusMotif = ''.join([max(item, key=item.count) for item in [*zip(*allBestMotifs)]])
        return allBestMotifs, consensusMotif, allBestScore

    def scrambleSequence(self):
        """
        scramble the DNA sequence to create a baseline
        Returns:
            a list of scrambled sequences

        """
        scrambleSequences = []
        for sequence in self.sequences:
            # print('o is', sequence)
            sequence = list(sequence)
            # https://www.w3schools.com/python/ref_random_shuffle.asp
            # shuffle is inplace
            shuffle(sequence)
            # print('s is', ''.join(sequence))
            scrambleSequences.append(''.join(sequence))
        # print(self.sequences)
        return scrambleSequences

    def iterScramble(self):
        """
        iterate randomized search many times with scrambled data
        Returns:
            the consensus motif and best score from scrambled data
        """
        original = self.sequences
        # print('original', original,'\n', self.sequences)
        self.sequences = self.scrambleSequence()
        _, consensusMotif, allBestScore = self.iterateSearch()
        # change the DNA sequences in the object back to original sequences
        self.sequences = original
        return consensusMotif, allBestScore

    def printResults(self):
        """
        Pretty printing the results
        """
        allBestMotifs, consensusMotif, allBestScore = self.iterateSearch()
        if self.scramble:
            scrambleConsensus, scrambleScore = self.iterScramble()
            print('Baseline Score with shuffled run result is', scrambleConsensus, scrambleScore)
            print('Non-shuffled run result is                ', consensusMotif, allBestScore)
        else:
            print('Non-shuffled run result is', consensusMotif, allBestScore)

        if self.matrix:
            for i in range(len(self.headers)):
                print(self.headers[i], allBestMotifs[i])

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
        # print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do
    except Usage as err:
        print(err.msg)

    # Get the commandline arguments
    iterations = myCommandLine.args.iterations
    motifLength = myCommandLine.args.motifLength
    pseudocount = myCommandLine.args.pseudocount

    scramble = myCommandLine.args.scramble
    gibbs = myCommandLine.args.gibbs
    matrix = myCommandLine.args.matrix

    # print('extracted args', iterations, motifLength, pseudocount, scramble, gibbs, matrix)

    fastaFile = FastAreader().readFasta()
    # store all sequence in a list
    headers = []
    sequences = []
    for header, sequence in fastaFile:
        headers.append(header.split()[0])
        sequences.append(sequence)
        # print(headers)
        # print(sequences)
    # print('DNA seqs are',sequences, 'length is', len(sequences))

    searchDNA = SearchMotif(headers, sequences, iterations, motifLength, pseudocount, scramble, gibbs, matrix)
    searchDNA.printResults()

if __name__ == "__main__":
    # start = time.time()
    main()
    # print('time consumed is', time.time() - start)
