#!/usr/bin/env python3
# Name: William Sobolewski
# Group Members: Jacob Curren, Dylan Brown, Zoe Petroff

from sequenceAnalysis import NucParams
from sequenceAnalysis import FastAreader
import sys

def main (fileName=None):
    """Calls FastAreader and NucParams, prints output"""
    myReader = FastAreader(fileName)# Create an instance of FastAreader with the given file name
    myNuc = NucParams() # Create an instance of NucParams
    for head, seq in myReader.readFasta() :# Iterate through the fasta sequences and add each sequence to myNuc
        myNuc.addSequence(seq)
    
    #print sequence length in Mb
    nucCount = (myNuc.nucCount()) / 1000000 #convert to Mb
    print("sequence length = {0:.2f} Mb".format(nucCount)) # 2 decimals
    
    #print GC content as a percentage
    GCcontent = ((myNuc.nucComp['G']+ myNuc.nucComp['C']) / (myNuc.nucCount())) * 100 #calculates percentage of G and C by getting value from nucComp dict
    print('\nGC content = {0:.1f}%\n '.format(GCcontent)) #prints GC content
     
    # Print the codon usage and frequency for each amino acid 
    for codon, aa in sorted(myNuc.rnaCodonTable.items(),key=lambda x: x[1]): # sorts the amino acids in alphabetical order
        if aa in myNuc.aaComp:
            aaCounts = myNuc.aaComp[aa]
            codonCounts = myNuc.codonComp[codon]
            val = (codonCounts/aaCounts)
            print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val*100, codonCounts))



if __name__ == "__main__":
    main() #calls main function with sys.stdin as file name "python genomeAnalyzer.py <filname.fa"