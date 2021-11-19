# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:55:48 2020

@author: Xavier Hernandez-Alias
"""

import os
import sys
import numpy as np
import pandas as pd

def translate_codon(codon):
    GENETIC_CODE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    codon = codon.replace("-","")
    if codon in GENETIC_CODE.keys():
        aa = GENETIC_CODE[codon]
        return aa
def split_codons(humanseq):
    n=0
    codpositions = [0]
    while len(humanseq) > n:
        until = 3
        while (until - humanseq[n:n+until].count("-"))<3:
            until += 1
        n = n+until
        codpositions.append(n)
    return codpositions


    
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
fasta_files = {file[:-3]:os.path.join(fasta_dir,batch,file) for batch in os.listdir(fasta_dir) if os.path.isdir(os.path.join(fasta_dir,batch)) for file in os.listdir(os.path.join(fasta_dir,batch)) if file.endswith(".fa")}

gc1df = pd.DataFrame()
gc2df = pd.DataFrame()
gc3df = pd.DataFrame()
for file in fasta_files.keys():
    # Get cds and gene info
    splitted_name = file.split("_")
    cds = splitted_name[1]
    gc1df.loc[cds,"gene"] = splitted_name[0]
    gc2df.loc[cds,"gene"] = splitted_name[0]
    gc3df.loc[cds,"gene"] = splitted_name[0]
    # Open file
    alignment = openfasta(fasta_files[file])
    # Compare each codon to that of human
    codpositions = split_codons(alignment["Human"])
    humancod = [alignment["Human"][codpositions[n]:codpositions[n+1]] for n in range(len(codpositions)-1)]
    humanaa = [translate_codon(c) for c in humancod]
    for species in alignment.keys():
        speciescod = [alignment[species][codpositions[n]:codpositions[n+1]] for n in range(len(codpositions)-1)]
        speciesaa = [translate_codon(c) for c in speciescod]
        gc1 = [cod[0] in ["G","C"] for n,cod in enumerate(speciescod) if (humanaa[n]==speciesaa[n]) and humanaa[n] and speciesaa[n]]
        gc2 = [cod[1] in ["G","C"] for n,cod in enumerate(speciescod) if (humanaa[n]==speciesaa[n]) and humanaa[n] and speciesaa[n]]
        gc3 = [cod[2] in ["G","C"] for n,cod in enumerate(speciescod) if (humanaa[n]==speciesaa[n]) and humanaa[n] and speciesaa[n]]
        if len(gc1)>=10:
            gc1df.loc[cds,species] = sum(gc1)/len(gc1)
            gc2df.loc[cds,species] = sum(gc2)/len(gc2)
            gc3df.loc[cds,species] = sum(gc3)/len(gc3)

# Save results
gc1df.to_csv(os.path.join("CDS_GC1.tsv"), sep="\t")
gc2df.to_csv(os.path.join("CDS_GC2.tsv"), sep="\t")
gc3df.to_csv(os.path.join("CDS_GC3.tsv"), sep="\t")
