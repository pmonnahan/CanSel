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

# Define Functions
def format_results(stat, stat_name, chrom, windows, pop):


    starts = np.take(windows, 0, axis=1)
    stops = np.take(windows, 1, axis=1)
    new_dat = pd.DataFrame(
        {'pop': [pop for x in range(0, len(stat))], 'chrom': [chrom for x in range(0, len(stat))],
         'start': starts, 'stop': stops, 'variable': [stat_name for x in range(0, len(stat))], 'value': stat})
    return new_dat


def get_ACdata(zarr_folder, chrom, samples, start=-9, stop=-9):
    callset = zarr.open_group(zarr_folder, mode='r')

    pos = callset[chrom]['variants']['POS']

    gt = allel.GenotypeDaskArray(callset[str(chrom)]['calldata']['GT'])  # Retrieve genotype data
    gt = gt.take(samples, axis=1).compute()  # subset data to samples of interest

    ac = gt.count_alleles()

    if start == -9:  start = min(pos)
    if stop == -9:  stop = max(pos)

    return ac, pos, start, stop

def get_biallelic(zarr_folder, chrom, samples):
    callset = zarr.open_group(zarr_folder, mode='r')

    gt = allel.GenotypeDaskArray(callset[str(chrom)]['calldata']['GT'])  # Retrieve genotype data
    gt = gt.take(samples, axis=1).compute()  # subset data to samples of interest

    ac = gt.count_alleles()
    sites = ac.is_biallelic_01()[:]

    return sites



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
    parser.add_argument('-o', type=str, metavar='output', required=True, help='')
    args = parser.parse_args()

    if args.c == 'all':
        chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']
    else:
        chroms = args.c.split(",")

    callset = zarr.open_group(args.z, mode='r')
    samples = callset[chroms[0]]['samples']
    winsize = int(args.w)

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

    if args.p2 != "None":
        try:
            loc2_samples = df[df.Population_code == args.p2].callset_index.values
        except IndexError:
            print("Did not find population " + args.p2 + " in popKey file.")
            args.p2 = "None"

    df_list = []

    for chrom in chroms:
        # pdb.set_trace()
        for pop in pops:
            loc_samples = df[df.Population_code == pop].callset_index.values
            ac, pos, start, stop = get_ACdata(args.z, chrom, loc_samples)
            if any(k in args.s for k in ['dxy', 'FD', 'Fst_Hud', 'Fst_Pat']) and args.p2 != "None":
                ac2, pos, start, stop = get_ACdata(args.z, chrom, loc2_samples)
                biallelic = get_biallelic(args.z, chrom, np.concatenate((loc_samples, loc2_samples), axis=0))
                ac = ac.compress(biallelic, axis=0)[:, :2]
                ac2 = ac2.compress(biallelic, axis=0)[:, :2]
                pos = pos.get_mask_selection(biallelic)
            if 'pi' in args.s:
                pi, windows, n_bases, counts = allel.windowed_diversity(
                    pos, ac, size=winsize, start=start, stop=stop, step=winsize / 2)
                new_dat = format_results(stat=pi, stat_name="pi", chrom=chrom, windows=windows, pop=pop)
                df_list.append(new_dat)
            if 'tajD' in args.s:
                tajD, windows, counts = allel.windowed_tajima_d(
                    pos, ac, size=winsize, start=start, stop=stop, step=winsize / 2)
                new_dat = format_results(stat=tajD, stat_name="tajD", chrom=chrom, windows=windows, pop=pop)
                df_list.append(new_dat)
            if 'dxy' in args.s and args.p2 != "None":
                dxy, windows, n_bases, counts = allel.windowed_divergence(
                    pos, ac, ac2, size=winsize, start=start, stop=stop, step=winsize / 2)
                new_dat = format_results(stat=dxy, stat_name="dxy", chrom=chrom, windows=windows, pop=pop)
                df_list.append(new_dat)
            if 'FD' in args.s and args.p2 != "None":
                FD, windows, n_bases, counts = allel.windowed_df(pos, ac, ac2, size=winsize, start=start, stop=stop,
                                                                 step=winsize / 2)
                new_dat = format_results(stat=FD, stat_name="FD", chrom=chrom, windows=windows, pop=pop)
                df_list.append(new_dat)
            if 'Fst_Hud' in args.s and args.p2 != "None":
                Fst_Hud, windows, counts = allel.windowed_hudson_fst(pos, ac, ac2, size=winsize, start=start,
                                                                     stop=stop, step=winsize / 2)
                new_dat = format_results(stat=Fst_Hud, stat_name="Fst_Hud", chrom=chrom, windows=windows, pop=pop)
                df_list.append(new_dat)
            if 'Fst_Pat' in args.s and args.p2 != "None":
                Fst_Pat, windows, counts = allel.windowed_patterson_fst(pos, ac, ac2, size=winsize, start=start,
                                                                        stop=stop, step=winsize / 2)
                new_dat = format_results(stat=Fst_Pat, stat_name="Fst_Pat", chrom=chrom, windows=windows, pop=pop)
                df_list.append(new_dat)

            if any(k in args.s for k in ['iHS', 'nSL']):
                gt, pos = get_hapdata(args.z, chrom, loc_samples)
                if args.P == 'true':
                    polarize = polarizeHapStat(args.z, chrom)
                else:
                    polarize = [1 for j in callset[chrom]['variants']['POS']]
            if 'iHS' in args.s:
                ihs, pos = calciHS(gt, pos, args.m)
                ihs = [x * j for x, j in zip(ihs, polarize)]
                new_dat = pd.DataFrame(
                    {'pop': [pop for x in range(0, len(ihs))], 'chrom': [chrom for x in range(0, len(ihs))],
                     'start': pos, 'stop': pos, 'variable': ['iHS' for x in range(0, len(ihs))], 'value': ihs})
                df_list.append(new_dat)
            if 'nSL' in args.s:
                nsl, pos = calcnSL(args.z, loc_samples, chrom)
                nsl = [x * j for x, j in zip(nsl, polarize)]
                new_dat = pd.DataFrame(
                    {'pop': [pop for x in range(0, len(nsl))], 'chrom': [chrom for x in range(0, len(nsl))],
                     'start': pos, 'stop': pos, 'variable': ['nsl' for x in range(0, len(nsl))], 'value': nsl})
                df_list.append(new_dat)

    if df_list:
        results = pd.concat(df_list)
        results.to_csv(args.o)
    else:
        print "no results to write"



