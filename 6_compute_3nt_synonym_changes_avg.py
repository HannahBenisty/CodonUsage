# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:55:48 2020

@author: Xavier Hernandez-Alias
"""

import os
import sys
import pandas as pd

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

def translate_codon(codon):
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
        
#%% Read fastas and compute codon identity

fasta_dir = sys.argv[1]

# Get list of fasta files
fasta_files = {file[:-3]:os.path.join(fasta_dir,batch,file) for batch in os.listdir(fasta_dir) if os.path.isdir(os.path.join(fasta_dir,batch)) for file in os.listdir(os.path.join(fasta_dir,batch)) if file.endswith(".fa")}
nfastas = len(fasta_files)

tuples_index = [(c1,c2) for c1 in GENETIC_CODE.keys() for c2 in GENETIC_CODE.keys() if (GENETIC_CODE[c1]==GENETIC_CODE[c2]) and (c1[0:2]==c2[0:2]) and (c1[2]!=c2[2])]
multi_index = pd.MultiIndex.from_tuples(tuples_index, names=['From', 'To'])
outdf = pd.DataFrame(index = multi_index)
for n,file in enumerate(fasta_files.keys()):
    if n%100==0:
        print("Progress: %i out of %i" % (n,nfastas))
    # Get cds and gene info
    splitted_name = file.split("_")
    cds = splitted_name[1]
    tempdf = pd.DataFrame(index = multi_index)
    # Open file
    alignment = openfasta(fasta_files[file])
    # Compare each codon to that of human
    codpositions = split_codons(alignment["Human"])
    humancod = [alignment["Human"][codpositions[n]:codpositions[n+1]] for n in range(len(codpositions)-1)]
    humanaa = [translate_codon(c) for c in humancod]
    for species in alignment.keys():
        if species not in tempdf.columns:
            tempdf[species] = 0
        speciescod = [alignment[species][codpositions[n]:codpositions[n+1]] for n in range(len(codpositions)-1)]
        speciesaa = [translate_codon(c) for c in speciescod]
        identity = [(cod.replace("-",""),humancod[n].replace("-","")) for n,cod in enumerate(speciescod) if (cod.replace("-",""),humancod[n].replace("-","")) in multi_index and humanaa[n] and speciesaa[n]]
        if len(identity)>=10:
            for c in identity:
                tempdf.loc[c,species] += 1
    # Compute percentage of changes by codon
    percdf = pd.DataFrame(index = multi_index, columns = tempdf.columns)
    cod1_in_df = list(set([s[0] for s in multi_index]))
    for c1 in cod1_in_df:
        cod2_in_df = [(s[0],s[1]) for s in multi_index if s[0]==c1]
        total = tempdf.loc[cod2_in_df,:].sum(axis=0)
        if total.mean()>=3:
            percdf.loc[cod2_in_df,:] = tempdf.loc[cod2_in_df,:].div(total)
    outdf[cds] = percdf.mean(axis=1)

# Save results
outdf.to_csv(os.path.join(fasta_dir,"CDS_3nt_synonym_changes_avg.tsv"), sep="\t")
