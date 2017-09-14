#!/usr/bin/env python3

import argparse

# this script will parse NGS index files and determine if index swapping occured
#looks at index reads determines if proper unique dual index or indexes were swapped and keeps a count per occurence
# filters quality scores on a per nucleotide basis if NT below cuttoff read is not retained
# reports undetermined indexes due to sequenceing error and such 

#requires 6 arguments to run listed below

# deifine argparse arguments
def get_arguments():
	parser = argparse.ArgumentParser(description="open  index read file(2) check q score of index reads and if above min cuttoff + valid index combo retain counts .")
	parser.add_argument("-f1", "--infile1", help="1st index read file (FASTQ format (FASTQ format", required=True, type=str)
	parser.add_argument("-f2", "--infile2", help="2nd index read file (FASTQ format", required=True, type=str)
	parser.add_argument("-q", "--qualityScoreCutoff", help="specifiy the quality score cutoff for index reads", required=True, type=int)
	parser.add_argument("-t", "--outfile", help="Name of tsv output(.tsv) file of index combos plus their counts", required=True, type=str)
        parser.add_argument("-u", "--undetermined_outfile", help="Name of tsv output(.tsv) file of undetermined index combos that pass qual filtering plus their counts", required=True, type=str)
	parser.add_argument("-f5", "--infile5", help="list of all known index tab seperated", required=True, type=str)
	return parser.parse_args()
args=get_arguments()

#intialize input file var
file1=args.infile1 
file2=args.infile2
file5=args.infile5

 #intialize output file var
fileOut=args.outfile
ud_fileOut=args.undetermined_outfile 

index={} #create a dict to hold all combos of index values and their counts
undetermined={} # create a dict to hold undetermined index combos that pass quality filtering and their counts

l=[] #create a list to hold index reads from index.tsv file ie this list holds all the possible index combinations that should be seen.

qual=args.qualityScoreCutoff # quality score cuttoff ie 30 

retained=False #boolean to indicate whether a record is kept or not (see below)

#open indexes.txt strip and split by tab grab the index reads and put in list l **this will be a list of all the indexes
with open(file5) as tsvfile:
    for line in tsvfile:
        parts=line.strip().split("\t")
        x=(parts[4])
        l.append(x)
del l[0] # the first line of this list is the header so remove that *** comment this line out if file has no header line 

#generate all possible combonations of 2 index values and store in dict with value @ 0
# ie index1,index2:0
for i in l:
    for j in l:
        index[i,j]=0

#function to convert phred score
def convert_phred(letter):
    """Converts a single character into a phred score"""
    
    n = ord(letter) - 33
    return n

#function to get the reverse complement of the index **paired end data only
def ReverseComplement1(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])


#open both index files and read and store the 4 indivdual lines to new variables as long as not EOF
with open(file1,"r")as R1,open(file2,"r")as R2:
    while True:
        
        header1=R1.readline()
        header2=R2.readline()
        if header1=="": #EOF  **note all files the same length so only need one EOF assert
            break
        seq1=R1.readline().strip()
        seq2=R2.readline().strip() 
        plus1=R1.readline().strip()
        plus2=R2.readline().strip()
        qual1=R1.readline().strip()
        qual2=R2.readline().strip()
        
        seq2=ReverseComplement1(seq2) # get the reverse comp of seq2 ** in order to proparly compare to seq1
        
        val=seq1,seq2 #variable to check if in index dict (see below)
    
        
        #read the quality score lines and convert the phred scores
        for i,c in zip(qual1,qual2): # **note zip is a built in to read multiple files simultainiously
            r=convert_phred(i) #ie convert qual 1 index read1
            g=convert_phred(c) #ie convert qual 2 index read2
            
            #if the quality score of each nucleotide in the index read is greater than min cutoff keep the record else dont
            if r >= qual and g >= qual:
                retained=True
            else:
                retained=False
                break
        #if keeping the record        
        if retained==True:
            if val in index: #if seq2,seq3 combo is in the combo dictionary
                index[val]+=1 #increase the count for that combo
                
                
            #if combo not in combo dict append the index to undetermined dict **ie seq error and such
            else:
                 if val in undetermined: #if seq1,2 val in undetermined dict increment count by 1
                    undetermined[val] += 1
                 else:
                    undetermined[val] = 1 #else put it in and count 1
                    
#output to tsv all index combos and their counts
with open(fileOut,"w") as f: #output as a tsv
        for key, value in index.items():
                f.write("{}\t{}\n".format(key, value))
                        
#output to tsv all undetermined index combos and their counts

with open(ud_fileOut,"w") as ud: #output as a tsv
        for key, value in undetermined.items():
                ud.write("{}\t{}\n".format(key, value))