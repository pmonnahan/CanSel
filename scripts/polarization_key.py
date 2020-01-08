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

# Define Functions
def polarizeHapStat(zarr_folder, chrom, capitalize=True, remove_missing_aa=True):
    # TODO: add polarization of iHS results by flipping sign when ref is derived. What do "-" and "." mean!?!?!
    # These take far too long to create de novo each time.  Should make once and
    callset = zarr.open_group(zarr_folder, mode='r')
    AA = callset[chrom]['variants']['AA']
    REF = callset[chrom]['variants']['REF']

    missing = 0
    unknown = 0
    edge = 0
    polarize = {}
    for j, aa in enumerate(AA):
        # if j % 1000 == 0: print(j)
        ref = REF[j]
        aa = aa.split("|")[0]
        alleles = ['A', 'G', 'C', 'T']
        id = callset[chrom]['variants']['ID'][j]
        if ref in alleles:  # establishes SNP variant
            alleles.remove(ref)
            if aa in alleles:  # ancestral allele differes from ref allele
                polarize[id] = -1
                print callset[chrom]['variants']['ID'][j]
            # if aa == ref or aa == ".":  # Reference allele matches annotated ancestral allele
            #     polarize.append(1)
            elif ref.lower() == aa:  # ancestral allele inference is low quality
                if capitalize:
                    polarize[id] = 1
                    # polarize.append(1)
                if not capitalize:
                    # polarize.append(0)
                    polarize[id] = 0
                    missing += 1
            # elif aa in alleles:  # Ancestral allele is found in list of non-ref allele possibilities
            #     polarize.append(-1)
            elif aa.capitalize() in alleles:
                if capitalize:
                    polarize[id] = -1
                    print callset[chrom]['variants']['ID'][j]
                else:
                    polarize[id] = 0
                    missing += 1
            elif 'unknown' in aa:
                if remove_missing_aa:
                    polarize[id] = 0
                    missing += 1
                else:
                    polarize[id] = 1
                unknown += 1
            else:
                polarize[id] = 0
                edge += 1
    print "Found " + str(unknown) + " variants with no known ancestral allele annotation."
    if remove_missing_aa:
        print "Converted unknowns to missing"
    print "Found " + str(missing) + " variants to be converted to missing."
    print "Found " + str(edge) + " edge cases that were set as missing"
    return(polarize)

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='')
    parser.add_argument('--use_lower_case', help='include lower case (less certain) ancestral allele calls ',
                        default=False, action="store_true")
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    args = parser.parse_args()

    for chrom in [x for x in range(1,23)] + ['X', 'Y']:
        pol_key = polarizeHapStat(args.i, chrom, capitalize=args.use_lower_case)
        pickle.dump(pol_key, open(args.o, 'wb'))

