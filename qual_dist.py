#!/usr/bin/env python3


# this script will parse a fastq file and determine the per line average quality score and the per nucleotide quality score 
# output is  average per nucleotide score for the entire file as a tsv and average per read qual score with counts as a tsv

# function"""Converts a single character into a phred score"""
def convert_phred(letter):
        n = ord(letter) - 33
        return n
i=0 # count for number of records calc
mean_scores = [0.0] *101  # Dict to hold the per NT Qual score **note change length based on read length
lineDict={} # dict to hold average line quality score and their counts

# open the file grab the quality score line, sum the line quality score and individual NT score
with open("R1_test.fastq","r") as fh:
   
    for line in fh:
        line=line.strip("\n") #strip new line char
        i+=1
        #if i%500000==0: #test to make sure its running
            #print(i/4)
        sum_lineQual=0
        if i%4==0: #grab the qual score line
            j=0
            for c in line: #go through each char in a line
                x=convert_phred(c) # convert score
                sum_lineQual +=x #keep a sum of scores for that line
                mean_scores[j] = mean_scores[j] +x # store per nucleotide score in mean scores ie each bp position running total
                j=j+1
            
            meanLineQual=sum_lineQual/len(line) #calculate the mean for the current line
            #print (meanLineQual)
            if meanLineQual in lineDict: #if mean line quality score in dict increment count by 1
                lineDict[meanLineQual] += 1
            else:
                lineDict[meanLineQual] = 1 #else put it in and count 1
        #i+=1
NR = i/4 #how many records

# convert all summed per bp qual scores to means for that bp postion and output to a tab sep file
with open("meanBasePairQual.tsv","w") as mq:
    
    for i in range(101): # **note change the range based on length of read 
        mean_scores[i]=mean_scores[i]/NR
        print(i,"\t",mean_scores[i],file=mq)

#export per line qual scores w counts to a tab sep file    
with open("meanLineQual.tsv","w") as f:
    for key, value in lineDict.items(): # ****This is how to write dictionary to tsv
        f.write("{}\t{}\n".format(key, value))
                
#print(mean_scores)
#print(NR)
