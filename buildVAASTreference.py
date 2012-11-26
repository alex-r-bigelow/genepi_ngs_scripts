#!/usr/bin/env python
import argparse, os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Leaves only chromosomes 1-22,X in a .vcf file')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .fasta file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .fa file')
    
    args = parser.parse_args()
    
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    for line in infile:
        if line.startswith(">"):
            chromosome = line[1:].split()[0]
            outfile.write(">chr%s\n"%chromosome)
        else:
            outfile.write(line)
    infile.close()
    outfile.close()