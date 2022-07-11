# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:55:48 2020

@author: user
"""

import pandas as pd
import numpy as np

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

#%% Read fasta and compute GC

# Get fasta
introns = openfasta("/home/xhernandez/Downloads/genome_assemblies_genome_fasta/ncbi-genomes-2022-04-19/GRCh38.p12.introns.fna")

# Build dataframe
GCdf = pd.DataFrame(index = pd.MultiIndex.from_tuples([], names=("ccds","gene","chr")),
					columns = ["#G","#C","#A","#T"])
for intron,seq in introns.items():
	ccds = intron.split(":")[0].split("_")[0]
	gene = intron.split(":")[0].split("_")[1]
	chrom = intron.split(":")[0].split("_")[2]
	if (ccds,gene,chrom) in GCdf.index:
		GCdf.loc[(ccds,gene,chrom),"#G"] += seq.upper().count("G")
		GCdf.loc[(ccds,gene,chrom),"#C"] += seq.upper().count("C")
		GCdf.loc[(ccds,gene,chrom),"#A"] += seq.upper().count("A")
		GCdf.loc[(ccds,gene,chrom),"#T"] += seq.upper().count("T")
	else:
		GCdf.loc[(ccds,gene,chrom)] = [seq.upper().count("G"),seq.upper().count("C"),
								 seq.upper().count("A"),seq.upper().count("T")]

# Compute percent
GCdf["GCcontent"] = np.divide(GCdf[["#G","#C"]].sum(axis=1),GCdf[["#G","#C","#A","#T"]].sum(axis=1))

#%% Save results
GCdf.to_csv("results/CCDSintrons_GCcontent.tsv", sep="\t")
