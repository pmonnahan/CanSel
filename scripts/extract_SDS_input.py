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
#os.environ["NUMEXPR_MAX_THREADS"]="272"
import pdb
import pandas as pd
import zarr
import numpy as np
import pickle

# Define Functions
def get_singletons(zarr_folder, chrom, samples, start=-9, stop=-9):
    callset = zarr.open_group(zarr_folder, mode='r')

    pos = callset[chrom]['variants']['POS']
    # pdb.set_trace()
    ref = callset[chrom]['variants']['REF']
    alt = callset[chrom]['variants']['ALT']
    ids = callset[chrom]['variants']['ID']

    gt = allel.GenotypeDaskArray(callset[str(chrom)]['calldata']['GT'])  # Retrieve genotype data
    gt = gt.take(samples, axis=1).compute()  # subset data to samples of interest

    ac = gt.count_alleles()

    if start == -9:  start = min(pos)
    if stop == -9:  stop = max(pos)

    flt = ac.is_singleton(1)
    pos2 = pos.get_mask_selection(flt)
    gf = gt.compress(flt, axis=0)
    sing_dict = {p:i for p,i in zip(pos2, np.where(gf.is_het())[1])}
    ind_dict = {}
    for key, value in sing_dict.items():
        if value in ind_dict:
            ind_dict[value].append(key)
        else:
            ind_dict[value] = [key]

    return ind_dict, gt, ids, ref, alt, pos, start, stop

#g = allel.GenotypeChunkedArray(callset[chrom]['calldata']['genotype'])
## the '[:]' syntax pulls the data from compressed storage into a numpy array
#ac = g.count_alleles()[:] #differs for zarr array
#Count singletons
#np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1))
#flt = ac.is_singleton(1)
#gf = g.compress(flt, axis=0)

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-z', type=str, metavar='zarr_file', required=True, help='')
    parser.add_argument('-k', type=str, metavar='popKey_file', required=True, help='file associating samples in VCF '
                                                                                   'to populations')
    parser.add_argument('-p', type=str, metavar='populations', default='all', help='comma-separated list (no spaces) '
                                                                                   'of populations that you want to '
                                                                                   'calculate diversity for.  Use '
                                                                                   '"all" to calculate for all '
                                                                                   'populations')
    parser.add_argument('-c', type=str, metavar='chromosomes', default='all', help='comma-separated list (no spaces) '
                                                                                   'of chromosomes that you want to '
                                                                                   'calculate diversity for.  Use '
                                                                                   '"all" to calculate for all '
                                                                                   'chromosomes')
    parser.add_argument('-P', type=str, metavar='polarization_key', default='None',
                        help='key for polarizing haplotype stat results; generated by polarization_key.py')
    parser.add_argument('-m', type=float, metavar='minimum_frequency', default=0.05,
                        help='minimum frequency to be included in test snps')
    parser.add_argument('-s', type=int, metavar='snps_per_file', default=999999999999999,
                        help='Controls the number of SNPs to write to a file.  If this value is greater than the total number of SNPs, a single file will be written')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    args = parser.parse_args()



# Prequisites: 	   Assumes input files are sorted by genomic coordinates and have a single matched chromosome.
# 		   Assumes the order of individuals matches for all input files.
# Arguments:	   s_file  	   The singletons.
# 		   		   A space/tab delimited file, that can be gzipped.
# 		   		   Entry {row i, column j} is the genomic position of the j'th singleton of individual i.
# 				   Singleton positions are assumed to be sorted.
# 				   Individuals are assumed the same order as in the test-SNPs file (t_file).
# 				   The singletons file is loaded into memory.
# 		   t_file  	   The test-SNPs.
# 		   		   A space/tab delimited file, that can be gzipped.
# 				   Each row describes a test SNP.
# 				   First 3 entries are: <snp-id> <ancestral-allele> <derived-allele> <position>.
# 				   The reminder are the individual genotypes (0/1/2) by the derived allele (i.e., 0=A/A;1=A/D;2=D/D).
# 				   Order of individuals should be as in the singletons file (s_file).

    snps_per_file = args.s

    if args.c == 'all':
        chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']
    else:
        chroms = args.c.split(",")

    callset = zarr.open_group(args.z, mode='r')
    samples = callset[chroms[0]]['samples']
    min_freq = args.m

    if args.P != "None":
        if not os.path.exists(args.P):
            print("Did not find polarization key at " + args.P)
        else:
            pol_key = pickle.load(open(args.P, 'rb'))

    df = pd.read_csv(args.k, sep="\t", header=0) # Read in population key file
    df = df[df['Sample_name'].isin(samples)] # Remove samples from pop key that are not in VCF data
    if not np.all(samples == df['Sample_name'].values): # Check that order of samples is same in VCF and pop key
        samples_list = list(samples)
        samples_callset_index = [samples_list.index(s) for s in df['Sample_name']]
        df['callset_index'] = samples_callset_index

    if args.p == 'all':
        pops = list(df.Population_code.unique())
    else:
        pops = args.p.split(",")

    df_list = []

    for chrom in chroms:
        for pop in pops:
            loc_samples = df[df.Population_code == pop].callset_index.values
            if len(loc_samples) == 0:
                print("No samples found in VCF for pop " + pop)
                break
            sing, gt, ids, a_allele, d_allele, pos, start, stop = get_singletons(args.z, chrom, loc_samples)

            if args.P != "None":
                polarize = [pol_key[x] for x in
                            ids[...]]  # Creates list of items used to polarize in same order as var_IDs
            else:
                polarize = [1 for x in ids[...]]

            gt = gt.to_n_alt()

            with open(f"{args.o}.bound", 'w') as bounds_file:
                bounds_file.write("1 9999999999999999\n")

            file_num = 1
            snps_file = open(f"{args.o}.{file_num}.snps", 'w')
            obs_file = open(f"{args.o}.{file_num}.obs", 'w')
            snp_count = 0
            for i, j in enumerate(polarize):
                freq = sum(gt[i]) / gt[i].shape[0]
                if j != 0 and min_freq < freq < (1 - min_freq):
                    snp_count += 1
                    if j == 1:
                        GT = [str(x) for x in gt[i].tolist()]
                        cmd = f"{ids[i]} {a_allele[i]} {d_allele[i][0]} {pos[i]} {' '.join(GT)}\n"
                    elif j == -1:
                        GT = [str(2 - x) for x in gt[i].tolist()]
                        cmd = f"{ids[i]} {d_allele[i][0]} {a_allele[i]} {pos[i]} {' '.join(GT)}\n"

                    if snp_count % snps_per_file == 0:
                        obs_file.write("\n")
                        snps_file.close()
                        obs_file.close()
                        file_num += 1
                        snps_file = open(f"{args.o}.{file_num}.snps", 'w')
                        obs_file = open(f"{args.o}.{file_num}.obs", 'w')

                    snps_file.write(cmd)
                    obs_file.write(f"1 ")

            obs_file.write("\n")
            snps_file.close()
            obs_file.close()

            with open(f"{args.o}.sings", 'w') as singletons_file:
                for ind in range(0,gt.shape[1]):
                    ising = sing[ind]
                    ising.sort()
                    ising = [str(x) for x in ising]
                    singletons_file.write(" ".join(ising) + "\n")