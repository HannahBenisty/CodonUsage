# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 12:04:52 2020

@author: user
"""

import pandas as pd

# Import mapping file
mapping = pd.read_csv("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt",
                      sep= "\t")

# Keep only CCDS with "Public" status, and remove duplicated CCDS
mapping = mapping.loc[mapping.ccds_status=="Public",:].drop_duplicates(["cds_locations"])
# Get list of genes
exonbed = pd.DataFrame(columns=["chr","start","end","name","score","strand"])
intronbed = pd.DataFrame(columns=["chr","start","end","name","score","strand"])
ce = 0; ci = 0
for idx in mapping.index:
	chrid = mapping.loc[idx,"#chromosome"]
	chrom = mapping.loc[idx,"nc_accession"]
	strand = mapping.loc[idx,"cds_strand"]
	ccds = mapping.loc[idx,"ccds_id"]
	gene = mapping.loc[idx,"gene"]
	splitEX = mapping.loc[idx,"cds_locations"][1:-1].split(", ")
	for n,ex in enumerate(splitEX):
		start = int(ex.split("-")[0])
		end = int(ex.split("-")[1])
		exonbed.loc[ce] = [chrom,start,end+1,"%s_%s_chr%s_%i" % (ccds,gene,chrid,n),".",strand]
		ce+=1
	for n in range(len(splitEX)-1):
		start = int(splitEX[n].split("-")[1])
		end = int(splitEX[n+1].split("-")[0])
		intronbed.loc[ci] = [chrom,start+1,end,"%s_%s_chr%s_%i" % (ccds,gene,chrid,n),".",strand]
		ci+=1

#%% Save
exonbed.to_csv("results/CCDSexons.bed", sep="\t", header=False, index=False)
intronbed.to_csv("results/CCDSintrons.bed", sep="\t", header=False, index=False)
