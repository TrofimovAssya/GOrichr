import sys, cPickle as pickle, numpy as np, os
from math import log,exp
from fisher import fisher_exact_nc, lnbico, lnfact, lngamm,func,ancestry
from go_prepper import prepper,slimmer,geneset_association,traceback,prepGOs


print sys.argv[1].split(".txt")[0]

"""  Script to perform GO term enrichment
Python rendition of the thresholded Fisher test
fisher test with threshold to avoid bias in big groups
originally by Sebastien Lemieux
optimized for use with pypy """


if not (os.path.exists("/Users/trofimov/Documents/Stage1415/GOrichr/Gset.p") and :
    prepGOs()
    
GOs = pickle.load(open("/Users/trofimov/Documents/Stage1415/GOrichr/Gset.p"))
go_dict = pickle.load(open("/Users/trofimov/Documents/Stage1415/GOrichr/go_dict.p"))
goslims = pickle.load(open("/Users/trofimov/Documents/Stage1415/GOrichr/goslim.p"))



# loading the gene set file
# gset will contain an array with a gene name per coordinate            
with open(sys.argv[1]) as f:
    for i in f.readlines()
        gset.append(i.rstrip())

"""
The following loop itterates through all GO terms, builds the 2 by 2 matrix and performs
Fisher test with threshold.
it retains only the GO terms that are statistically significant (p-val < 0.05)
into a auto-sorted array.
"""

results = []
rstats = []
golist = set()

for i in range (0,len (GOs)):
    k = GOs.keys()[i]
    both = set(GOs[k]) & set(gset)
    """
    n11 = number of genes present in gset and GOs
    n21 = number of genes not in gset but in GOs
    n12 = number of genes in gset but not in GOs
    n22 = number of genes not in GOs and not in gset (21024-gset)
    """
    n11 = len(both)
    n12 = len(gset) - n11
    n21 = len(GOs[k]) - n11
    if len(GOs[k]) > 21024: continue
    n22 = 21024 - (n11 + n21 + n12)
    
    p = fisher_exact_nc(n11,n12,n21,n22,2)
    n = lookup(k,go_dict)
    n = n.split("\n")[0]
    if p<0.05:
        results.append([p,n,k,(str(n)+str(k)),both])
        golist.update(k)
    

results.sort()

n = len(results)
stat=True #status see later
    
x = sys.argv[1]

deletes = set()
for j in golist:
    ancestry(j,False,go_dict,goslims,deletes)
    
for i in range(0,len(results)):
    if results[i][2] not in deletes: 
        if stat==True: #we want to create a file only if there is an enrichment!
            f = open(str("GOterms_"+x),"w")
            stat=False
        f.write(str(results[i])+"\n")
f.close()
