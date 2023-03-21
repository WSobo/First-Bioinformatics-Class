#!/usr/bin/env python3
# Name: William Sobolewski
# Group Members: Jacob Curren, Dylan Brown, Zoe Petroff


'''This is a module containing 3 classes: NucParams, ProtParams, and FastaAreader.
genomeAnalyzer.py imports this module'''



class NucParams:
    """
    A class that computes the nucleotide, codon and amino acid compositions of a given DNA or RNA sequence.

    Attributes:
    rnaCodonTable (dict): A dictionary that maps RNA codons to their corresponding amino acids.
    dnaCodonTable (dict): A dictionary that maps DNA codons to their corresponding amino acids.

    codonCount (dict): A dictionary that stores the count of each codon in the sequence.
    sequence (dict): A dictionary that stores the codons of the input sequence in RNA form.
    nucComp (dict): A dictionary that stores the count of each nucleotide in the sequence.
    codonComp (dict): A dictionary that stores the count of each codon in the sequence.
    aaComp (dict): A dictionary that stores the count of each amino acid in the sequence.

    Methods:
    addSequence (inSeq): Takes an input sequence, formats it and stores the codons in the `sequence` attribute.
    aaComposition (): Computes the amino acid composition of the sequence.
    nucComposition (): Computes the nucleotide composition of the sequence.
    codonComposition (): Computes the codon composition of the sequence.
    nucCount (): Computes the total count of nucleotides in the sequence.
    """
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C', # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-', # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W', # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    
    def __init__ (self, inString=''):
        self.codonCount = {}
        self.sequence = {}
        self.nucComp = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0}
        self.codonComp = {codon:0 for codon in NucParams.rnaCodonTable}
        self.aaComp = {AA:0 for AA in NucParams.rnaCodonTable.values()}
        if inString:
            self.addSequence(inString)
        
    def addSequence (self, inSeq):
        '''format the string and store the codons in self.sequence in RNA form'''
        for base in inSeq : #calculate composition of nucleotide bases
            if base in 'ACGTUN':
                self.nucComp[base] += 1
        
        for i in range(0, len(inSeq), 3): #for each codon in the sequence, convert to rna codons and determine the composition of codons and amino acids
            codon = inSeq[i:i+3] 
            codon = codon.replace('T', 'U')
            if codon in self.codonComp: #add one to the key value in self.codonComp and the associated aa to the codon in self.aaComp
                self.codonComp[codon] += 1
                self.aaComp[self.rnaCodonTable[codon]] += 1
    

    def aaComposition(self):
        """
        Compute the amino acid composition of the sequence

        Returns:
            A dictionary where the keys are the single letter amino acid codes and the values are the frequency of each amino acid in the sequence
        """
        return self.aaComp

    def nucComposition(self):
        """
        Compute the nucleotide composition of the sequence {ATGCNU}
        """
        return self.nucComp

    def codonComposition(self):
        """Computes the codon composition of the sequence."""
        return self.codonComp

    def nucCount(self):
        """Computes the total count of nucleotides in the sequence."""
        return sum(self.nucComp.values())

import collections
import bisect

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        self.protein = protein #initializes the input protein to be referenced in the definitions within the class
        

    def aaCount (self):
        '''returns the length of the amino acid string (not all that useful tbh)'''
        return len(self.protein)

    def pI (self, precision=2):
        """"This will return a pH that the input protein will be at a net neutral charge"""
        #for pH in range(0,1401): #upto pH of 14. 2 decimal places
            #charge = self._charge_(pH/100) #calls charge method
            #if (abs(charge) < 0.01): #2 decimal of pH
              #  return pH/100 #returns pH with 2 decimals
    
        #return -1 #returns -1 if no pH works maybe false would be better???
        
        
        '''This is the preforms the same as above, but uses the bisect package. EXTRA CREDIT'''
        low = 0   #low pH range value
        high = 14 #high pH range value
        
        while high - low > 10 ** -precision: #preforms a binary search using bisect. precision = 2 so to two decimal places
            mid = (low + high) / 2
            charge = self._charge_(mid) #charge at the average of the two pH values
            if charge > 0: #plays a game of higher lower
                low = mid
            elif charge < 0:
                high = mid
            else:
                return mid #once the optimal pH is determined (charge=0), return the value
        return (low + high) / 2

    def aaComposition (self) :
        '''returns a .Counter() object of all 20 AA, and if the AA arent in self.protein, there will still be a key and value of 0'''
        amino_acid_count = collections.Counter(dict.fromkeys(ProteinParam.aa2mw, 0)) #create a counter object with all 20 amino acids with value of 0. utilize the aa2mw dictionary with all the AA keys
        c = collections.Counter(self.protein) # Count how many of each AA are 
        return amino_acid_count + c #updates the values of all 
    
    def _charge_(self, pH):
        '''Does some crazy math for calculating the net charge of the AA string at specific pH values. This defintion is called by def pI()'''
        netCharge = 0 #initialize the object netCharge
        
        '''I split up the math for calculating the Positive and negative charges.'''
        for aa in self.aa2chargePos: #for each 'K', 'R', 'H', in the input protein 
            if aa in self.aaComposition() and aa in self.aa2chargePos: 
                netCharge += self.aaComposition()[aa] * (10 ** self.aa2chargePos[aa]) / (10 ** self.aa2chargePos[aa] + 10 ** pH) #do the math that was shown in the instructions provided and add to the netCharge
        netCharge += (10 ** self.aaNterm) / ((10 ** self.aaNterm ) + (10 ** pH)) # add on the charge of the single N-terminus

        for aa in self.aa2chargeNeg: #Now do the same for the Negative charged AA's
            if aa in self.aaComposition() and aa in self.aa2chargeNeg: #for each 'D', 'E', 'C', and 'Y' in the input protein
                netCharge -= self.aaComposition()[aa] * (10 ** pH) / (10 ** self.aa2chargeNeg[aa] + 10 ** pH) # do math and subtract from the netCharge

        netCharge -= (10 ** pH) / ((10 ** self.aaCterm) + (10 ** pH)) # subtract the charge of the single C-terminus

        return netCharge #return the netCharge (utilized by the definition pI())


    def molarExtinction (self):
        '''This definition utilizes dark magic to determine a coefficient that can be utilized for determining the concentration of the protein when suspended in water'''
        mE = 0 #intializes molar extinction coefficent
       
        '''This magic spell involves utilizing the fact that the amino acids: 'Y', 'W', and 'C' all absorb light (to varying degrees) at 280nm. these next lines are just counting how many of each of the three are in the protein'''
        Ycount = 0
        Wcount = 0
        Ccount = 0
        for i in range(len(self.protein)):
            if self.protein[i] in ProteinParam.aa2abs280:
                if self.protein[i] == 'Y':
                    Ycount += 1
                elif self.protein[i] == 'W':
                    Wcount += 1
                elif self.protein[i] == 'C':
                    Ccount += 1
                    
        '''Now that the amount of each have been counted, they're multiplied by their asosiated absorption value listed in the dictionary aa2abs280, and added together'''
        mE = (Ycount * ProteinParam.aa2abs280['Y']) + (Wcount * ProteinParam.aa2abs280['W']) + (Ccount * ProteinParam.aa2abs280['C'])
    
        return mE #returns the molar extinction coefficent

    def massExtinction (self):
        '''returns the mass extinction coefficent by dividing the molar extinction coefficient by the molecular weight of the AA string'''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0 #to prevent an error for dividing by zero incase the AA input is empty

    def molecularWeight (self):
        '''Returns the molecular weight of the protein by adding together the mass of each individual acid, then subtracting the water molecules that are removed during the peptide backbone formation'''
        weight = 0
        for i in range(len(self.protein)):
            weight += ProteinParam.aa2mw[self.protein[i]] #references the aa2mw dictionary provided to access the mass of each AA
        return weight - ((len(self.protein)-1) * self.mwH2O) #returns the molecularWeight after subtracting the mass of the # of waters 
        
    
import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen(self):
        """function to handle either a file name or a file-like object"""
        if self.fname == None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        print('checking self.fname:', self.fname)
        """
        try:
            fileH = self.doOpen()
            print('File handle:', fileH)
        except Exception as e:
            print('Error opening file:', e)
            return
        """
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            #starts reading the DNA sequence
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
