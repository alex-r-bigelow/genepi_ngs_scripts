#!/usr/bin/env python
import argparse
from genome_utils import chromosomesToKeep

def run(args):
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    for line in infile:
        if line.startswith("#"):
            if 'contig' in line:
                contigID = line[line.find('ID=')+3:]
                before = line[:line.find('ID=')+3]
                after = contigID[contigID.find(','):]
                contigID = contigID[:contigID.find(',')]
                
                if not contigID.startswith('chr'):
                    contigID = 'chr' + contigID
                if contigID in chromosomesToKeep:
                    outfile.write(before + contigID + after)
            else:
                outfile.write(line)
                continue
        columns = line.strip()
        if len(columns) <= 1:
            continue
        columns = columns.split('\t')
        chrom = columns[0]
        if not chrom.startswith('chr'):
            chrom = "chr" + chrom
        if chrom in chromosomesToKeep:
            outfile.write("%s\t%s\n"%(chrom,'\t'.join(columns[1:])))
    infile.close()
    outfile.close()

if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description='Leaves only chromosomes 1-22,X,Y in a .vcf file')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    
    args = parser.parse_args()
    run(args)