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
import pandas as pd
import zarr
import numpy as np
import pickle


def get_hapdata(zarr_folder, chrom, samples, min_maf=0.05):
    callset = zarr.open_group(zarr_folder, mode='r')
    pos = callset[chrom]['variants']['POS']
    gt = allel.GenotypeDaskArray(callset[str(chrom)]['calldata']['GT'])  # Retrieve genotype data
    gt = gt.take(samples, axis=1).compute()  # subset data to samples of interest

    # Filter for allele frequency
    ac = gt.count_alleles()
    af = ac.to_frequencies()
    af = np.take(af, 0, axis=1)
    keep = [True if (1 - min_maf > x > min_maf) else False for x in af]

    gt = gt.compress(keep, axis=0)
    pos = [x for x in pos[:] * keep if x != 0]

    return gt, pos

def format_results(stat, stat_name, chrom, pos, pop, polarization_key, ac="None", n_bins=20):

    stat = [x * j for x, j in zip(stat, polarization_key)]
    if ac != "None":
        allel.standardize_by_allele_count(stat, ac, n_bins=n_bins)
        stat_name += ".std"
    new_dat = pd.DataFrame(
        {'pop': [pop for x in range(0, len(stat))], 'chrom': [chrom for x in range(0, len(stat))],
         'start': pos, 'stop': pos, 'variable': [stat_name for x in range(0, len(stat))], 'value': stat})
    return new_dat


def calciHS(gt, pos, map_file, map_cols=['chrom', 'bp', 'cM'], map_sep=r"\s+"):

    if map_file == "None":
        haps = gt.to_haplotypes()
        ihs = allel.ihs(haps, pos)
    else:
        map = pd.read_csv(map_file, delimiter=map_sep, names=map_cols, header=None, low_memory=False)
        map = map[map.chrom == chrom]
        # min_map = max(map.bp[map.cM==min(map.cM)])  # Remove positions that are outside of bp range in gen_map
        # max_map = min(map.bp[map.cM==max(map.cM)])  # Need to account for ties in cM at distal positions
        #
        # keep = [True if (min_map < x < max_map) else False for x in pos]
        # gt = gt.compress(keep, axis=0)
        # pos = [x for x in np.array(pos) * np.array(keep) if x != 0]
        new_map = np.interp(pos, map.bp, map.cM)
        haps = gt.to_haplotypes()
        ihs = allel.ihs(haps, pos, map_pos=new_map)

    return ihs


def calcXPEHH(gt, gt2, pos, map_file, map_cols=['chrom', 'bp', 'cM'], map_sep=r"\s+"):

    if map_file == "None":
        h1 = gt.to_haplotypes()
        h2 = gt2.to_haplotypes()
        xpehh = allel.xpehh(h1, h2, pos, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
    else:
        map = pd.read_csv(map_file, delimiter=map_sep, names=map_cols, header=None, low_memory=False)
        map = map[map.chrom == chrom]
        # min_map = max(map.bp[map.cM==min(map.cM)])  # Remove positions that are outside of bp range in gen_map
        # max_map = min(map.bp[map.cM==max(map.cM)])  # Need to account for ties in cM at distal positions
        #
        # keep = [True if (min_map < x < max_map) else False for x in pos]
        # gt = gt.compress(keep, axis=0)
        # pos = [x for x in np.array(pos) * np.array(keep) if x != 0]
        new_map = np.interp(pos, map.bp, map.cM)
        h1 = gt.to_haplotypes()
        h2 = gt2.to_haplotypes()
        xpehh = allel.xpehh(h1, h2, pos, map_pos=new_map, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)

    return xpehh


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-z', type=str, metavar='zarr_file', required=True, help='')
    parser.add_argument('-k', type=str, metavar='popKey_file', required=True, help='file associating samples in VCF '
                                                                                   'to populations')
    parser.add_argument('-p', type=str, metavar='populations', default='all', help='comma-separated list (no spaces) '
                                                                                   'of populations that you want to '
                                                                                   'calculate stats for.  Use '
                                                                                   '"all" to calculate for all '
                                                                                   'populations')
    parser.add_argument('-p2', type=str, metavar='comparison_population', default='CEU', help='which population to '
                                                                                              'use for comparison')
    parser.add_argument('-c', type=str, metavar='chromosomes', default='all', help='comma-separated list (no spaces) '
                                                                                   'of chromosomes that you want to '
                                                                                   'calculate diversity for.  Use '
                                                                                   '"all" to calculate for all '
                                                                                   'chromosomes')
    parser.add_argument('-s', type=str, metavar='statistics', default='all', help='comma-separated list (no spaces) '
                                                                                  'of statistics to be calculated.  '
                                                                                  'Options are [pi, tajD, iHS, nSL]')
    parser.add_argument('-w', type=str, metavar='window_size', default='50000', help='step size is 1/2 winsize')
    parser.add_argument('-m', type=str, metavar='genetic_map_file', default='None', help='')
    parser.add_argument('-P', type=str, metavar='polarization_key', default='None', help='key for polarizing haplotype stat results; generated by polarization_key.py')
    parser.add_argument('--standardize', help='standardize haplotype scores within allele frequency bins ',
                        default=False, action="store_true")
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    args = parser.parse_args()

    if args.p == 'all':
        pops = list(df.Population_code.unique())
    else:
        pops = args.p.split(",")

    if args.c == 'all':
        chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']
    else:
        chroms = args.c.split(",")

    callset = zarr.open_group(args.z, mode='r')
    samples = callset[chroms[0]]['samples']
    winsize = int(args.w)

    # TODO: integrate and debug output of polarization_key.py for use in this script
    # TODO: add haplotype diversity calculations under name 'hapDiv'. no windowed approach available in allel currently
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

    if args.p2 != "None":
        try:
            loc2_samples = df[df.Population_code == args.p2].callset_index.values
        except IndexError:
            print("Did not find population " + args.p2 + " in popKey file.")
            args.p2 = "None"

    df_list = []

    for chrom in chroms:
        pdb.set_trace()
        if any(k in args.s for k in ['xpehh', 'xpnSL']) and args.p2 != "None":
            gt2, pos = get_hapdata(args.z, chrom, loc2_samples)
        var_IDs = callset[chrom]['variants/ID']
        if args.P != "None":
            polarize = [pol_key[x] for x in var_IDs]
        else:
            polarize = [1 for x in var_IDs]
        for pop in pops:
            loc_samples = df[df.Population_code == pop].callset_index.values
            gt, pos = get_hapdata(args.z, chrom, loc_samples)
            if args.standardize:
                ac = gt.count_alleles()
                if args.p2 != "None":
                    ac2 = gt2.count_alleles()
                    ac = np.add(ac, ac2)
            else:
                ac = "None"
            if 'iHS' in args.s:
                ihs = calciHS(gt, pos, args.m)
                new_dat = format_results(ihs, 'iHS', chrom, pos, pop, polarize, ac)
                df_list.append(new_dat)
            if 'nSL' in args.s:
                nsl = allel.nsl(gt.to_haplotypes(), pos)
                new_dat = format_results(nsl, 'nSL', chrom, pos, pop, polarize, ac)
                df_list.append(new_dat)
            if 'xpehh' in args.s and args.p2 != "None":
                xpehh = calcXPEHH(gt, gt2, pos, args.m)
                new_dat = format_results(xpehh, 'xpehh', chrom, pos, pop, polarize, ac)
                df_list.append(new_dat)
            if 'xpnSL' in args.s and args.p2 != "None":
                xpnsl = allel.xpnsl(gt.to_haplotypes(), gt2.to_haplotypes())
                new_dat = format_results(xpnsl, 'xpnSL', chrom, pos, pop, polarize, ac)
                df_list.append(new_dat)

    if df_list:
        results = pd.concat(df_list)
        results.to_csv(args.o)
    else:
        print "no results to write"