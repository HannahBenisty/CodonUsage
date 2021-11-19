# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 16:18:50 2020

@author: Xavier Hernandez-Alias
"""

import os
import sys
    
def openfasta(filename):
    # Open the file
    fileHandle = open(filename)
    # Read all the lines in the file
    seqs = fileHandle.read().split(">") # Split file whenever a > is detected
    # Close the file
    fileHandle.close()
    # Loop over the lines
    seqdict = {}
    for seq in seqs:
        if seq: # Skip empty sequences
            splitted = seq.split("\n") # Split different lines of one same sequence
            # the first line is the name of the sequence, the rest is the sequence itself
            seqdict[splitted[0]] = "".join(splitted[1:])
    return seqdict

#%% Read fastas and compute codon identity

fasta_dir = sys.argv[1]

# Get list of fasta files
fasta_files = {file[:-3]:os.path.join(fasta_dir,batch,file) for batch in os.listdir(fasta_dir) if os.path.isdir(os.path.join(fasta_dir,batch)) for file in os.listdir(os.path.join(fasta_dir,batch)) if file.endswith(".fa")}
humanfasta = open(os.path.join("human_CDS.fa"),"w") 
for file in fasta_files:
    # Open file
    alignment = openfasta(fasta_files[file])
    # Compare each codon to that of human
    humanseq = alignment["Human"].replace("-","")
    # Write sequence
    humanfasta.write(str(">%s\n") % file)
    humanfasta.write(str("%s\n") % humanseq)

humanfasta.close()
