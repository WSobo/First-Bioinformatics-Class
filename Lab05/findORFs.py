#to use this code, run in command line from the directory this file and sequenceAnalysis.py and the fasta file are in

#cd C:\Users\Wcsob\OneDrive\Desktop\BME\BME160\Lab05

#then and example of what you could input and parameters for output:

#python findORFs.py <tass2.fa >tass2ORFdata-ATG-100.txt -mG 300 -lG -s ATG -t GTG



import sys
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

def main(inFile = None, options = None):
    from sequenceAnalysis import FastAreader, orfFinder
    '''
    Find some genes.  
    '''
    #parse the command line options
    thisCommandLine = CommandLine(options)
    reader = FastAreader(inFile)#Create a FastAreader object to read the input file
    sequences = reader.readFasta()
    
    ###### replace the code between comments.
    print (thisCommandLine.args)
    # thisCommandLine.args.longestGene is True if only the longest Gene is desired
    longestGene = thisCommandLine.args.longestGene
    # thisCommandLine.args.start is a list of start codons
    start = thisCommandLine.args.start
    # thisCommandLine.args.stop is a list of stop codons
    stop = thisCommandLine.args.stop
    # thisCommandLine.args.minGene is the minimum Gene length to include
    minGene = thisCommandLine.args.minGene
    
    
    #######
    
    for header, sequence in sequences:
        of = orfFinder(sequence, header, minGene, longestGene, start, stop,)
        of.findORFs() #finds all the orfs on the top strand and appends to pGlist
        of.compStrand() #finds all the orfs on the bottom strand and appends to pGlist
        of.outputFile() #sorts by length and adds all the identified orfs to the specified output file
    
if __name__ == "__main__":
    main() #runs the code