#!/bin/python/
import sys
import cPickle as pickle
import numpy as np
from math import log,exp
import os.path
import os



        
### fisher test with threshold to avoid bias in big groups
#### originally by Sebastien Lemieux
### optimized for use with pypy


def lngamm(z):
    ## Reference: "Lanczos, C. 'A precision approximation 
    ## of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
    ## Translation of  Alan Miller's FORTRAN-implementation
    ## See http://lib.stat.cmu.edu/apstat/245

    x = 0.1659470187408462e-06 / (z + 7);
    x += 0.9934937113930748e-05 / (z + 6);
    x -= 0.1385710331296526     / (z + 5);
    x += 12.50734324009056      / (z + 4);
    x -= 176.6150291498386      / (z + 3);
    x += 771.3234287757674      / (z + 2);
    x -= 1259.139216722289      / (z + 1);
    x += 676.5203681218835      / (z);
    x += 0.9999999999995183;
    return (log(x) - 5.58106146679532777 - z + (z - 0.5) * log (z + 6.5) )


def lnfact (n):
    if n<=1:
        return (0)
    return (lngamm(n+1))
    

def lnbico (n,k):
    return (lnfact (n) - lnfact (k) - lnfact (n-k))

def fisher_exact_nc (n11,n12,n21,n22,w):
    ##Fisher's exact test modified to use Fisher's non-central hypergeometric
    ##distribution with a odds-ratio bias of w.  The procedure returns the p-value
    ##based on a null-hypothesis of the odds-ratio being <= w.
    ##Significant calls indicates that n11 / n12 is enriched by at least w.
    x = n11
    m1 = n11 + n21
    m2 = n12 + n22
    n = n11 + n12
    x_min = max(0, n - m2)
    x_max = min(n, m1)
    l = map(lambda y: (lnbico (m1, y) + lnbico (m2, n-y) + y* log (w)),range(x_min,x_max+1))
    max_l = max(l)
    sum_l = log(reduce (lambda x,y: x+y, map(lambda x: exp(x-max_l), l),0))
    den_sum = log(reduce(lambda x,y: x+y, map(lambda x: exp(l [ x-x_min] - max_l), range(x,x_max+1)) ,0))
    
    return exp (den_sum - sum_l)

def func(k,go_dict):
    g = go_dict[k]
    return g["name"]

### function ancestry, makes sure to keep the highest enriched ancestor
def ancestry (child,slim,go_dict,goslims,deletes):
    if go_dict.has_key(child):
        cnames = go_dict[child]
        for i in cnames.keys():
            if i[:4] == "is_a":
                parent = cnames[i].split(" ! ")[0]
                if (parent in golist):
                    if (slim == False):
                        deletes.update([child])
                    if (goslims.has_key(parent)):
                        ancestry(parent,True,go_dict,goslims,deletes)
                        
                    else:
                        ancestry(parent,False,go_dict,goslims,deletes)

    


