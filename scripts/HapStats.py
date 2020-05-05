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
from math import isnan
import scipy

# TODO: COMMENT!!

def get_hapdata(zarr_folder, chrom, samples, biallelic_sites):
    callset = zarr.open_group(zarr_folder, mode='r')
    pos = callset[chrom]['variants']['POS']
    gt = allel.GenotypeDaskArray(callset[str(chrom)]['calldata']['GT'])  # Retrieve genotype data
    gt = gt.take(samples, axis=1).compute()  # subset data to samples of interest
    gt = gt.compress(biallelic_sites, axis=0)
    pos = pos.get_mask_selection(biallelic_sites)

    return gt, pos


def get_biallelic(zarr_folder, chrom, samples, min_maf = 0.01):
    callset = zarr.open_group(zarr_folder, mode='r')

    gt = allel.GenotypeDaskArray(callset[str(chrom)]['calldata']['GT'])  # Retrieve genotype data
    gt = gt.take(samples, axis=1).compute()  # subset data to samples of interest

    ac = gt.count_alleles()
    af = ac.to_frequencies()
    af = np.take(af, 0, axis=1)
    minMAF = [True if (1 - min_maf > x > min_maf) else False for x in af]
    # pdb.set_trace()
    biallelic = ac.is_biallelic_01()[:]
    sites = [x * y for x, y in zip(minMAF, biallelic)]

    return sites

def format_results(stat, stat_name, chrom, pos, pop, polarization_key, ac="None", n_bins=20):

    stat = [x * j for x, j in zip(stat.tolist(), polarization_key) if j != 0]  # Polarization will set stat values to 0 if ancestral allele is unknown
    pos = [x for x, j in zip(pos.tolist(), polarization_key) if j != 0]  # These 0 sites need to be removed from all arrays because they crash standardization
    if ac != "None":
        flt = np.array([True if x != 0 else False for x in polarization_key])
        ac = ac.compress(flt, axis=0)
        stat, bins = allel.standardize_by_allele_count(stat, ac[:, 1], n_bins=n_bins, diagnostics=False)
        stat_name += ".std"
    pos = [y for x, y in zip(stat.tolist(), pos) if not isnan(x)]
    stat = [x for x in stat.tolist() if not isnan(x)]
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
        new_map = np.interp(pos, map.bp, map.cM)
        h1 = gt.to_haplotypes()
        h2 = gt2.to_haplotypes()
        xpehh = allel.xpehh(h1, h2, pos, map_pos=new_map, min_ehh=0.05, include_edges=True, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)

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
    parser.add_argument('-p2', type=str, metavar='comparison_population', default='None', help='which population to '
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

    if args.s == 'all':
        args.s = ['iHS', 'nSL', 'xpehh', 'xpnSL']
        if args.p2 == "None":
            args.s = [x for x in args.s if x not in ['xpehh', 'xpnSL']]
    else:
        args.s = args.s.split(",")
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

    df = pd.read_csv(args.k, sep="\t", header=0)  # Read in population key file
    df = df[df['Sample_name'].isin(samples)]  # Remove samples from pop key that are not in VCF data
    if not np.all(samples == df['Sample_name'].values):  # Check that order of samples is same in VCF and pop key
        samples_list = list(samples)
        samples_callset_index = [samples_list.index(s) for s in df['Sample_name']]
        df['callset_index'] = samples_callset_index
    if args.p == 'all':
        pops = list(df.Population_code.unique())
    else:
        pops = args.p.split(",")

    if args.p2 != "None":
        try:
            loc2_samples = df[df.Population_code == args.p2].callset_index.values
            if len(loc2_samples) == 0:
                print("No samples found in VCF for reference pop " + args.p2)
                args.p2 = "None"
        except IndexError:
            print("Did not find population " + args.p2 + " in popKey file.")
            args.p2 = "None"


    df_list = []
    # pdb.set_trace()
    for chrom in chroms:
        # pdb.set_trace()

        for pop in pops:
            loc_samples = df[df.Population_code == pop].callset_index.values
            if len(loc_samples) == 0:
                print("No samples found in VCF for pop " + pop)
                break
            var_IDs = callset[chrom]['variants/ID']
            if any(k in args.s for k in ['xpehh', 'xpnSL']) and args.p2 != "None":
                if pop == args.p2:
                    print("Must specify two different populations for between population metrics")
                    break
                else:
                    biallelic = get_biallelic(args.z, chrom, np.concatenate((loc_samples, loc2_samples), axis=0))
                    gt, pos = get_hapdata(args.z, chrom, loc_samples, biallelic)
                    gt2, pos = get_hapdata(args.z, chrom, loc2_samples, biallelic)
                    var_IDs = var_IDs.get_mask_selection(biallelic)
            else:
                biallelic = get_biallelic(args.z, chrom, loc_samples)
                gt, pos = get_hapdata(args.z, chrom, loc_samples, biallelic)
                var_IDs = var_IDs.get_mask_selection(biallelic)

            if args.P != "None":
                polarize = [pol_key[x] for x in
                            var_IDs[...]]  # Creates list of items used to polarize in same order as var_IDs
            else:
                polarize = [1 for x in var_IDs[...]]

            if args.standardize:
                ac = gt.count_alleles()
                if args.p2 != "None":
                    ac2 = gt2.count_alleles()
                    ac = np.add(ac, ac2)  # does not work
            else:
                ac = "None"

            if 'iHS' in args.s:
                ihs = calciHS(gt, pos, args.m)
                new_dat = format_results(ihs, 'iHS', chrom, pos, pop, polarize, ac)
                df_list.append(new_dat)
            if 'nSL' in args.s:
                nsl = allel.nsl(gt.to_haplotypes())
                new_dat = format_results(nsl, 'nSL', chrom, pos, pop, polarize, ac)
                df_list.append(new_dat)
            if 'xpehh' in args.s:
                if args.p2 == "None":
                    print("Must specify a second population via -p2 in order to calculate xpehh")
                else:
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
        print("no results to write")
        oo = open(args.o, 'w')  # Create empty output file for snakemake bookkeeping
        oo.close()
