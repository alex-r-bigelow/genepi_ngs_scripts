#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Applies "DP Filter" to variants that are "5 or 6 sigma from the '+
                                 'mean coverage across all samples" (see GATK best practice v4')
parser.add_argument('--sigma', type=int, dest="sigma",
                    help='number of standard deviations below which a variant will PASS')
parser.add_argument('--in', type=str, dest="in",
                    help='input .vcf file')
parser.add_argument('--out', type=str, dest="out",
                    help='output .vcf file')

args = parser.parse_args()
print args.sigma
print args.vcf