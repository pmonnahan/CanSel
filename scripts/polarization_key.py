#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: 
Requirements: 
 Takes the following arguments:
    -i: 
    -o:
Date:
"""

# Import Modules
import allel
import subprocess
import argparse
import os
import pdb
import zarr
import pickle
import line_profiler

# Define Functions
# @profile
def polarizeHapStat(zarr_folder, chrom, capitalize=True, remove_missing_aa=True):
    # TODO: add polarization of iHS results by flipping sign when ref is derived. What do "-" and "." mean!?!?!
    # These take far too long to create de novo each time.  Should make once and
    callset = zarr.open_group(zarr_folder, mode='r')
    AA = callset[chrom]['variants']['AA'][...]
    REF = callset[chrom]['variants']['REF'][...]
    IDS = callset[chrom]['variants']['ID'][...]

    counts = [0,0,0]
    polarize = {}
    for j, aa in enumerate(AA):
        if j % 100 == 0: print(j)
        ref = REF[j]
        aa = aa.split("|")[0]
        alleles = ['A', 'G', 'C', 'T']
        id = IDS[j]
        # pdb.set_trace()
        if aa == '.' or aa == '':
            polarize[id] = 0
            counts[0] += 1
        elif ref in alleles:  # establishes SNP variant
            alleles.remove(ref)
            if ref == aa:  # ancestral allele matches reference
                polarize[id] = 1
                counts[2] += 1
            elif ref.lower() == aa:  # ancestral allele inference is low quality but matches reference
                if capitalize:
                    polarize[id] = 1
                    counts[2] += 1
                if not capitalize:
                    polarize[id] = 0
                    counts[0] += 1
            elif aa in alleles:  # ancestral allele differes from ref allele
                polarize[id] = -1
                counts[1] += 1

            elif aa.capitalize() in alleles:
                if capitalize:
                    polarize[id] = -1
                    counts[1] += 1
                else:
                    polarize[id] = 0
                    counts[0] += 1
            elif 'unknown' in aa:
                if remove_missing_aa:
                    polarize[id] = 0
                    counts[0] += 1
                else:
                    polarize[id] = 1
            else:
                polarize[id] = 0
                counts[0] += 1
    print "Found " + str(counts[0]) + " variants to be converted to missing."
    print "Found " + str(counts[1]) + " variants to be polarized"
    print "Found " + str(counts[2]) + " variants to be unchanged"
    return(polarize)

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='')
    parser.add_argument('--use_lower_case', help='include lower case (less certain) ancestral allele calls ',
                        default=False, action="store_true")
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    parser.add_argument('-c', type=str, metavar='chromosome', default="22", help='')
    args = parser.parse_args()

    if args.c == 'all':
        for chrom in [x for x in range(1,23)] + ['X', 'Y']:
            pol_key = polarizeHapStat(args.i, chrom, capitalize=args.use_lower_case)
            pickle.dump(pol_key, open(args.o, 'wb'))
    else:
        pol_key = polarizeHapStat(args.i, args.c, capitalize=args.use_lower_case)
        pickle.dump(pol_key, open(args.o, 'wb'))
