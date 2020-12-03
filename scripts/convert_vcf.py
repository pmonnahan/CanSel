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
import shutil
os.environ["NUMEXPR_MAX_THREADS"]="272"

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-vcf', type=str, metavar='vcf_file', required=True, help='')
    parser.add_argument('-name', type=str, metavar='vcf_name', required=True, help='')
    parser.add_argument('-o', type=str, metavar='output_path', required=True, help='Output_path')
    args = parser.parse_args()

    outdir = os.path.dirname(args.o)

    if not os.path.exists(outdir): os.mkdir(outdir)

    print("compressing chromosome " + args.name)
    allel.vcf_to_zarr(args.vcf, args.o, group=str(args.name), fields='*', overwrite=True)
    print("done compressing chromosome " + args.name)
