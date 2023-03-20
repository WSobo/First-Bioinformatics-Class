import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
#from Bio import pairwise2

"""
finds SNPS in a FASTA file of a Yersinia plasmid inputted via sys.stdin



Usage: finding mutations

Command line Example: python BME160Final.py <mutPlasmid.fa
"""



class yersiniaSNP:
    def __init__(self, refFile):
        self.referenceSeq = SeqIO.read(refFile, "fasta")
        self.querySeq = SeqIO.read(sys.stdin, "fasta")
        self.snpPositions = []
        
    def findMutation(self):
        """
        Compare the reference and query sequences base by base and find any differences.

        For each base pair in the sequences, compare the corresponding bases in the
        reference and query sequences. If there is a difference, record the position
        of the difference in the query sequence and add it to the list of SNP positions.
        Prints to sys.stdout

        Returns:
        None
        """
        for i in range(min(len(self.referenceSeq), len(self.querySeq))):
            refBase = self.referenceSeq[i] # Get the base from the reference sequence at index i
            queryBase = self.querySeq[i] # Get the base from the query sequence at index i
            if refBase != queryBase: # If the bases are different, we have a SNP
                refPos = i + 1 # Calculate the reference position of the SNP
                queryPos = (i + 1) // 70 + 2  # Calculate the line number of the query sequence where the SNP occurred
                # Print the SNP position and add it to the list of SNP positions
                print(f"SNP at line {queryPos} of the query FASTA file at base position: {refPos}: {refBase} -> {queryBase}")
                self.snpPositions.append(queryPos) #list is updated with the query position of the SNP.

                
    def displayGraph(self):
        """
        Display a graph showing the distribution of SNPs per line in the query FASTA file.
        If there are no SNPs, print a message to indicate that no mutations were found.
        """
        if not self.snpPositions:
            print("No mutations found.")
            return
        num_snps = [0] * 953
        for pos in self.snpPositions:
            num_snps[pos-1] += 1
        plt.subplots(figsize=(48, 6))
        plt.plot(range(953), num_snps)
        plt.xlabel("Line number in query FASTA file")
        plt.ylabel("Number of SNPs")
        plt.ylim(-1, max(num_snps)+1)
        plt.xlim(0, 960)
        plt.xticks(range(0, 953, 25))
        plt.title("SNP Distribution per line in FASTA file")
        plt.tight_layout()
        plt.subplots_adjust(left=0.06, bottom=0.08, right=0.96, top=0.92, wspace=None, hspace=None)
        plt.show()


    def orfIdentifier(self):
        """
        possible feature in the future. It's a joke rn
        """
        answer = input("Do you want to know the ORFs in which the SNPs occured? [Y/N]").replace(" ", "").upper()
        if answer != "Y":
            return
        elif answer == "Y":
            print("Too bad lol")





# Example usage: compare two FASTA files
def main():
    finder = yersiniaSNP("refPlasmid.fa")
    finder.findMutation()
    finder.displayGraph()
    #finder.orfIdentifier()
    
if __name__ == "__main__":
    main()