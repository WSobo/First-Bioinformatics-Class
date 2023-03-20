# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# Name: William Sobolewski (1846541)
# Group Members: Zoe Petroff, Dylan Brown, Jacob Curren
"""
command line example: python findUnique.py < bos-tRNA.fa > output.fa
"""
import sys 

class tRna:
    '''
    1.reads a fasta file containing the sequences of all the tRNA's.
    2.creates a powerset of all subsequences for each tRNA sequence.
    3.identifies the unique subsequences in all the set's
    4.
    '''
    tRNAObjects = [] # holds each TRna object that has its own header, sequence and power set

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.powerSet = self.getPowerSet(sequence)
        self.uniqueSubs = set()
        self.essentialSubs = set()
        tRna.tRNAObjects.append(self)

    def getPowerSet(self, sequence):
        '''
        a powerset is a set of all the possible substrings of each length of a string
        example: "ABCD" powerset is "ABCD", "ABC", "BCD", "AB", "BC", "CD", "A", "B", "C", "D"
        '''
        powerSet = set()
        for i in range(len(sequence)):
            for j in range(i+1, len(sequence)+1):
                powerSet.add(sequence[i:j])
        return powerSet

    def buildUniques(self):
        ''' 
        Result: self.uniqueSubs will contain a sets of only the unique substrings 
        '''
        # Create an empty set to store all the subsets from other RNA objects
        otherSubsets = set()
        # Loop through all the RNA objects except the current one
        for eachRNA in tRna.tRNAObjects:
            if eachRNA is not self:
                # Add all the subsets from each RNA object to the set
                otherSubsets = otherSubsets.union(eachRNA.powerSet)
        # Find the difference between the current RNA object's subsets and the others'
        self.uniqueSubs = self.powerSet - otherSubsets

    def findEssentials(self):
        '''
        compares subsets to eachother and keeps only the essential ones to each set stored in self.essentialSubs
        '''
        
        self.nonEssential = set()# Create an empty set to store the non-essential subsets
        
        for subset in self.uniqueSubs:# Loop through all the unique subsets
            # Find the index and length of the subset in the sequence
            index = self.sequence.find(subset)
            length = len(subset)
            # Find the left and right boundaries of the subset
            left = index - 1
            right = index + length + 1
            # If the left boundary is valid, add the extended subset to the non-essential set
            if not left < 0:
                self.nonEssential.add(self.sequence[left:index+length])
            # If the right boundary is valid, add the extended subset to the non-essential set
            if not right < 0:
                self.nonEssential.add(self.sequence[index:right])
        
        self.essentialSubs = self.uniqueSubs - self.nonEssential # Find the difference between the unique subsets and the non-essential subsets


class FastAreader:
    ''' 
    Define objects to read FastA files.
    '''
    def __init__(self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen(self):
        """function to handle either a file name or a file-like object"""
        if self.fname == None:
            return sys.stdin
        else:
            return open(self.fname, encoding='utf-8')
        
    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            for line in fileH:
                if line.startswith('>'):
                    if header:  # check if header is not empty
                        yield header, sequence
                        header = ''
                        sequence = ''
                    header = line[1:].rstrip()
                else:
                    sequence += line.strip()
            if header:  # return the last sequence
                yield header, sequence

def main():
    tRnaFinder = FastAreader() # create FastAreader object that takes from sys.stdin

    tRnaObjects = [] # empty list that will store tRNA objects

    for head, seq in tRnaFinder.readFasta():   # Loop through each tRNA sequence in the input and create a tRNA object for each one
        tRnaObj = tRna(head.replace(' ', ''), seq.replace('_', '').replace('-', '').replace('.', '')) # makes powerSet for each trna sequence
        tRnaObjects.append(tRnaObj)

    sortTrnas = sorted(tRnaObjects, key=lambda t: t.header.replace(' ', '')) # sorts the tRNAs so they can be iterated through alphabetically by their header

    for eachtRNA in tRnaObjects: # Find the unique and essential subsets for each tRNA object
        eachtRNA.buildUniques()
        eachtRNA.findEssentials()

    # Write to stdout in UTF-8 encoding
    with open(sys.stdout.fileno(), mode='w', encoding='utf-8', closefd=False) as stdout:
        for everyTrna in sortTrnas:
            # Write the header and sequence of the tRNA
            stdout.write(everyTrna.header.replace(' ', '') + '\n')
            stdout.write(everyTrna.sequence + '\n')
            # Sort the essential subsets by their position in the sequence and write them out
            sortedtRna = sorted(everyTrna.essentialSubs, key=lambda x: everyTrna.sequence.find(x))
            for e in sortedtRna:
                position = everyTrna.sequence.find(e)
                stdout.write('{}{}\n'.format('.' * position, e))

if __name__=="__main__":
    main()