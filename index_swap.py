#!/usr/bin/env python3

import argparse
import os
import gzip

# deifine argparse arguments
def get_arguments():
	parser = argparse.ArgumentParser(description="open paired end fastq files(2) and index read file(2) check q score of index reads and if above min cuttoff + valid index combo retain record of all 4 files in 4 new fastqs sorted by index read .")
	parser.add_argument("-f1", "--infile1", help="1st paired end Sequence file (FASTQ format", required=True, type=str)
	parser.add_argument("-f2", "--infile2", help="1st index read file (FASTQ format", required=True, type=str)
	parser.add_argument("-f3", "--infile3", help="2nd index read file (FASTQ format", required=True, type=str)
	parser.add_argument("-f4", "--infile4", help="2nd paired end sequence file (FASTQ format", required=True, type=str)
	parser.add_argument("-q", "--qualityScoreCutoff", help="specifiy the quality score cutoff for index reads", required=True, type=int)
	parser.add_argument("-t", "--outfile", help="Name of tsv output(.tsv) file of index combos plus their counts", required=True, type=str)
	parser.add_argument("-f5", "--infile5", help="list of all known index tab seperated", required=True, type=str)
	parser.add_argument("-o1", "--outFolder1", help="name of FOLDER were you want to store correct dual indexed reads(file names generate automaticaly", required=True, type=str)
	parser.add_argument("-o2", "--outFolder2", help="name of FOLDER you want to store index swaped reads(file names generate automaticaly", required=True, type=str)
	parser.add_argument("-o3", "--outFolder3", help="name of FOLDER you want to store undetermined reads(file names generate automaticaly", required=True, type=str)
	parser.add_argument("-z", "--zipped", help="are the sequenced reads in seperate zipped folders", required=True, type=bool)
	return parser.parse_args()
args=get_arguments()


file1=args.infile1 #intialize input file var
file2=args.infile2
file3=args.infile3
file4=args.infile4
file5=args.infile5
fileOut=args.outfile #intialize output file var
#folder1=args.outFolder1
#folder2=args.outFolder2
#folder3=args.outFolder3
zipped=args.zipped # is sequence files zipped (bool)

#create output folders
fold1 = args.outFolder1+"/"
os.makedirs(os.path.dirname(fold1), exist_ok=True)

fold2 = args.outFolder2+"/"
os.makedirs(os.path.dirname(fold2), exist_ok=True)

fold3 = args.outFolder3+"/"
os.makedirs(os.path.dirname(fold3), exist_ok=True)

index={} #create a dict to hold all combos of index values
l=[] #create a list to hold index reads from index.tsv file

op=open # for opening the sequenced reads(see below) using to determine open vs gzip.open

#if seq files are zipped set op to gzip.open for opening seq reads
if zipped ==True:
    op=gzip.open

#open index.tsv strip and split by tab grab the index reads and put in list l
with open(file5) as tsvfile:
    for line in tsvfile:
        parts=line.strip().split("\t")
        x=(parts[1])
        l.append(x)
    #print(l)

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


qual=args.qualityScoreCutoff # quality score cuttoff
retained=False #boolean to indicate whether a record is kept or not (see below)

#open all 4 files and read and store the 4 indivdual lines to new variables as long as not EOF
with op(file1,"rt")as R1,op(file2,"rt")as R2,op(file3,"rt")as R3,op(file4,"rt")as R4:
    
    while True:
        
        header1=R1.readline()
        header2=R2.readline()
        header3=R3.readline()
        header4=R4.readline()
        if header1=="": #EOF  **note all files the same length so only need one EOF assert
            break
        seq1=R1.readline().strip()
        seq2=R2.readline().strip() 
        seq3=R3.readline().strip()
        seq4=R4.readline().strip() 
        plus1=R1.readline().strip()
        plus2=R2.readline().strip()
        plus3=R3.readline().strip()
        plus4=R4.readline().strip()
        qual1=R1.readline().strip()
        qual2=R2.readline().strip()
        qual3=R3.readline().strip()
        qual4=R4.readline().strip()
        val=seq2,seq3 #variable to check if in index dict (see below)
        head="_"+seq2+"_"+seq3 #variable to change header sequence (see below)
        head1=seq2+"_"+seq3 #variable to write output file (see below)
        
        #read the quality score lines and convert the phred scores
        for i,c in zip(qual2,qual3): # **note zip is a built in to read multiple files simultainiously
            r=convert_phred(i) #ie convert qual 2 index read1
            g=convert_phred(c) #ie convert qual 3 index read2
            
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
                
                #split the headers on space and take the left side (ie remove 2 N... keep @ header) append the index sequences to the header( ie head variable)
                header1,header2,header3,header4=header1.split(" ")[0]+head,header2.split(" ")[0]+head,header3.split(" ")[0]+head,header4.split(" ")[0]+head
                
                if seq2==seq3: #if indexes are identical(correct dual index) place in specified output folder
                    # open a file in correct output format and append everything to the file
                    #** note files are opened as append so don't run twice or you will double file size
                    with open(fold1+head1+"_R1.fastq","a")as rof1,open(fold1+head1+"_R2.fastq","a")as rof2,open(fold1+head1+"_R3.fastq","a")as rof3,open(fold1+head1+"_R4.fastq","a")as rof4:
                        print(header1,"\n"+seq1,"\n"+plus1,"\n"+qual1,file=rof1)
                        print(header2,"\n"+seq2,"\n"+plus2,"\n"+qual2,file=rof2)
                        print(header3,"\n"+seq3,"\n"+plus3,"\n"+qual3,file=rof3)
                        print(header4,"\n"+seq4,"\n"+plus4,"\n"+qual4,file=rof4)
                else:
                    with open(fold2+head1+"_R1.fastq","a")as rof1,open(fold2+head1+"_R2.fastq","a")as rof2,open(fold2+head1+"_R3.fastq","a")as rof3,open(fold2+head1+"_R4.fastq","a")as rof4:
                        print(header1,"\n"+seq1,"\n"+plus1,"\n"+qual1,file=rof1)
                        print(header2,"\n"+seq2,"\n"+plus2,"\n"+qual2,file=rof2)
                        print(header3,"\n"+seq3,"\n"+plus3,"\n"+qual3,file=rof3)
                        print(header4,"\n"+seq4,"\n"+plus4,"\n"+qual4,file=rof4)
            #if combo not in combo dict append the index seq to header and output sequences to undetermined outfile
            else:
                header1,header4=header1.split(" ")[0]+head,header4.split(" ")[0]+head
                with open(fold3+"Undetermined_index_pair_R1.fastq","a")as urof1,open(fold3+"Undetermined_index_pair_R4.fastq","a")as urof4:
                    print(header1,"\n"+seq1,"\n"+plus1,"\n"+qual1,file=urof1)
                    print(header4,"\n"+seq4,"\n"+plus4,"\n"+qual4,file=urof4)
#output to tsv all index combos and their counts
with open(fileOut,"w") as f: #output as a tsv
        for key, value in index.items():
                f.write("{}\t{}\n".format(key, value))
                
                
