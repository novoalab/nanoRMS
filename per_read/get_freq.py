#!/usr/bin/env python3
desc="""Return modification frequency calucated using features extracted by get_features.py

TBD:
- 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerów/Barcelona, 30/11/2020
"""

import os, sys, numpy as np, pandas as pd
from datetime import datetime
from get_features import VERSION, logger
from common_functions import load_data_reps, get_freq, KMeans, KNeighborsClassifier
import pickle

def get_regions(bed):
    """Load candidate positions from BED file"""
    regions = []
    for l in open(bed):
        ldata = l[:-1].split(",")
        chrom, s, e = ldata[:3]
        s, e = int(s), int(e)
        if len(ldata)<6: strand = "+"
        else: strand = ldata[5]
        for pos in range(s, e):
            regions.append((chrom, pos+1, strand))
    return regions

def get_freq_diff(outfn, fasta, control, sample, bed, minCov=10, nn=1, 
                  features=["si", "tr"], names=("kmeans", "knn"),
                  clfs=(KMeans(n_clusters=2), KNeighborsClassifier())):
    """Report modification frequency for positive sites from BED file"""
    # load regions
    regions = get_regions(bed)
    logger("Processing %s region(s)...\n"%len(regions))
    
    strains_unique = ["control", "sample"]
    strain2bam = {strains_unique[0]: control, strains_unique[1]: sample}
    strains = [*[strains_unique[0]]*len(control), *[strains_unique[1]]*len(sample)]
    
    # load data
    logger("Loading features...\n")
    region2data = load_data_reps(fasta, control+sample, regions, features, strains, strains_unique, nn=nn)

    pickle.dump(region2data, open(outfn, 'wb'))
    return 
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-1", "--control", nargs="+", help="control BAM file(s)")
    parser.add_argument("-2", "--sample", nargs="+", help="sample BAM file(s)")
    parser.add_argument("-o", "--output", default="mod_freq.tsv.gz", help="output name [%(default)s]")
    parser.add_argument("-b", "--bed", required=1, help="positions to test")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-m", "--mincov", default=10, type=int, help="reference FASTA file")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))
    # check parameters
    if o.mincov<10:
        sys.stderr.write("[ERROR] -m/--mincov must be at least 5!")
        sys.exit(1)
        
    # encode tombo output into BAM files
    get_freq_diff(o.output, o.fasta, o.control, o.sample, o.bed, o.mincov, 
                  nn=3)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s    \n"%dt)

