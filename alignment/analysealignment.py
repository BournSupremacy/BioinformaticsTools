#!/usr/bin/env python

def simplematches(alignmentfile, seqcount):
    
    # This function counts the number of matches that the query file had against the reference/target.
    # It does this by checking the 10th column of the psl file or the 13th column of the show-coords file, which is the query sequence name.
    # There are no limits on coverage and identity.
    # It returns the percentage of contigs from the query that matched with the reference.

    hits = []
    
    if ".psl" in alignmentfile: 
        for line in open (alignmentfile, 'r'):
            splitline = line.split('\t')
            hits.append(splitline[9].strip()) # adding the query contig ID to the list
    
    elif ".show-coords" in alignmentfile: 
        with open(alignmentfile, 'r') as f:
            for i in range(0,4):
                next(f) # skipping first three lines, which are the headers
            for line in f:
                splitline = line.split('\t')
                hits.append(splitline[12].strip()) # adding the query contig ID to the list
    
    hits.sort()
    uniquematches = len(set(hits)) # finding only unique ones
    percentagematch = (uniquematches/seqcount)*100
    
    print (f"Number of query contigs that matched with the reference: {uniquematches}")
    print (f"Percentage of query contigs that matched with the reference: {percentagematch:.3f}")

def matches_with_threshold(alignmentfile, seqcount, identitylimit, coveragelimit):
    
    # This function finds the number of matches but with a threshold for both identity and coverage.
    # It calculates the identity and coverage for each hit, and then only keeps the ones that are above the threshold.
    # Identity is the percentage of residues that match up in the local alignment (which is each line),
    # i.e., identity = total bp matches divided by the local alignment length.
    # Coverage is the percentage of the query sequence length that is included in the alignment,
    # i.e., coverage = total bp matched and unmatched divided by total query sequence length.
    # For the .psl file:
        # Identity = column 1 + column 3/column 1 + column 2 + column 3 (total matches/total alignment length)
        # Coverage = column 1 + column 2 + column 3/column 11 (total alignment length/query length)
    # For the .show-coords file:
        # Identity = column 7
        # Coverage = column 11
    # It also adds up the total bp in the local alignments and keeps track of how many sequences from the query had a match.

    seqsandlengths = {} # dict to keep track of sequences with hits and their alignment with the most bp
    
    if ".psl" in alignmentfile:

        identitylimit = identitylimit/100
        coveragelimit = coveragelimit/100
        
        for line in open (alignmentfile, 'r'):
            identity = 0 # reset the variables for each line
            coverage = 0

            splitline = line.split('\t')
            matches = float(splitline[0].strip()) + float(splitline[2].strip())
            alignmentlength = float(splitline[0].strip()) + float(splitline[1].strip()) + float(splitline[2].strip()) 
            querylength = float(splitline[10].strip())

            identity = matches/alignmentlength
            coverage = alignmentlength/querylength

            if identity >= identitylimit and coverage >= coveragelimit:
                seqid = splitline[9].strip() # query contig ID
                if seqid in seqsandlengths.keys(): # check whether sequence has alerady got an alignment
                    if matches > seqsandlengths.get(seqid): # if yes, check whether this alignment is longer
                        seqsandlengths[seqid] = matches # if yes, update alignment
                else:
                    seqsandlengths[seqid] = matches # add it to the dictionary
            
        identitylimit = identitylimit*100
        coveragelimit = coveragelimit*100
                    
    elif ".show-coords" in alignmentfile:
        with open(alignmentfile, 'r') as f:
            for i in range(0,4):
                next(f) # skipping first three lines, which are the headers
            for line in f:
                identity = 0 # reset the variables for each line
                coverage = 0
                
                splitline = line.split('\t')
                identity = float(splitline[6].strip())
                coverage = float(splitline[10].strip())
                alignmentlength = float(splitline[5].strip())
                matches = alignmentlength*(identity/100) # NUCmer doesn't give exact matches, so work it out based on alignment length and identity
                
                if identity >= identitylimit and coverage >= coveragelimit:
                    seqid = splitline[12].strip() # query contig ID
                    if seqid in seqsandlengths.keys(): # check whether sequence has alerady got an alignment
                        if matches > seqsandlengths.get(seqid): # if yes, check whether this alignment is longer
                            seqsandlengths[seqid] = matches # if yes, update alignment
                    else:
                        seqsandlengths[seqid] = matches # add it to the dictionary
    
    uniquematches = len(seqsandlengths)
    percentagematch = (uniquematches/seqcount)*100
    
    # Add up total number of bp that had exact matches in the alignments
    totalbp = 0
    for val in seqsandlengths.values():
        totalbp = totalbp + val
    
    totalbp = round(totalbp)

    # print (f"Number of query contigs that matched with the target at {identitylimit}% identity and {coveragelimit}% coverage: {uniquematches}")
    # print (f"Percentage of query contigs that matched with the target at {identitylimit}% identity and {coveragelimit}% coverage: {percentagematch:.2f}")
    # print (f"Total exact bp matches: {totalbp} bp = {(totalbp/1000000):.2f} Mbp")
    
    print (f"{identitylimit}\t{coveragelimit}\t{uniquematches}\t{percentagematch:.2f}\t{totalbp}\t{(totalbp/1000000):.2f}")


# A script to analyse the .psl file from an alignment using pblat or the .show-coords file from an alignment using nucmer
# The first argument given is the ENTIRE PATH of the psl or show-coords file of interest.
# The second argument is the ENTIRE PATH of the original query file - this will allow us to find the sequence count.
# The third and fourth arguments are percentage identity and percentage coverage respectively.
   
import sys

alignmentfile = sys.argv[1]
queryfile = sys.argv[2]

# Count the number of sequences in the query file
file = open(queryfile, "r")
data = file.read()
seqcount = data.count(">")

if len(sys.argv) == 3:
    simplematches (alignmentfile, seqcount)

else:
    identitylimit = float(sys.argv[3])
    coveragelimit = float(sys.argv[4])

    if identitylimit and coveragelimit < 1: # making sure the percentage is in percentage format; if they are already, just leave it
        identitylimit = identitylimit*100
        coveragelimit = coveragelimit*100

    matches_with_threshold (alignmentfile, seqcount, identitylimit, coveragelimit)
