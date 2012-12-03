#!/usr/bin/env python
import argparse

if __name__ == '__main__':
    chromosomesToKeep = set(['chr1',
                             'chr2',
                             'chr3',
                             'chr4',
                             'chr5',
                             'chr6',
                             'chr7',
                             'chr8',
                             'chr9',
                             'chr10',
                             'chr11',
                             'chr12',
                             'chr13',
                             'chr14',
                             'chr15',
                             'chr16',
                             'chr17',
                             'chr18',
                             'chr19',
                             'chr20',
                             'chr21',
                             'chr22',
                             'chrX',
                             'chrY'])
    
    parser = argparse.ArgumentParser(description='Leaves only chromosomes 1-22,X,Y in a .fasta file, formatted in a way that VAAST won\'t complain about')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .fasta file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .fa file')
    
    args = parser.parse_args()
    
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    writeLines = True
    for line in infile:
        if line.startswith(">"):
            chromosome = line[1:].split()[0]
            if not line.lower().startswith("chr"):
                chromosome = "chr" + chromosome
            if not chromosome in chromosomesToKeep:
                writeLines = False
            else:
                writeLines = True
                outfile.write(">%s\n"%chromosome)
        elif writeLines:
            outfile.write(line)
    infile.close()
    outfile.close()