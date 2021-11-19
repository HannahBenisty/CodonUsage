# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 12:04:52 2020

@author: Xavier Hernandez-Alias
"""

import sys
import os
import pandas as pd
import requests

# Import mapping file
mapping = pd.read_csv("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt",
                      sep= "\t", index_col="ccds_id")

# Keep only CCDS with "Public" status, and remove duplicated CCDS
mapping = mapping.loc[mapping.ccds_status=="Public",:].drop_duplicates(["cds_locations"])
# Get list of genes
genefile = sys.argv[1]
os.mkdir("alignments")
filesperbatch = 2000
n=0
with open(genefile) as infile:
    for line in infile:
        gene = line.rstrip()
        # Get gene features
        if gene in mapping.gene.values:
            for cds in mapping.loc[mapping.gene==gene,:].index:
                if n%filesperbatch==0:
                    tempdir = str("batch%i" % (n/filesperbatch))
                    os.mkdir(os.path.join("alignments",tempdir))
                cdsinfo = mapping.loc[cds,:]
                chrom = cdsinfo["#chromosome"]
                strand = cdsinfo["cds_strand"]
                coords = None
                for i in cdsinfo["cds_locations"][1:-1].split(", "):
                    start = int(i.split("-")[0])+1
                    end = int(i.split("-")[1])+1
                    if coords:
                        coords = coords + str("+chr%s:%i-%i" % (chrom,start,end))
                    else:
                        coords = str("chr%s:%i-%i" % (chrom,start,end))
                    
                # Get alignment. More ifno here: https://data.broadinstitute.org/compbio1/cav.php?controlsState=show&intervals=chr1%3A64830953-64831021&alnset=hg38_58&spliceSites=human&Help=
                url = str("https://data.broadinstitute.org/compbio1/cav.php?intervals={}&strand={}&maxCodons=100000&spliceSites=human&alnset=hg38&fastaOut=".format(coords, strand))
                myfile = requests.get(url)
                # Write file
                open(os.path.join("alignments",tempdir,str("%s_%s.fa" % (gene,cds))), 'wb').write(myfile.content)
                n+=1
        else:
            print("Ignoring gene %s: coding sequence is not annotated in CCDS." % gene)
