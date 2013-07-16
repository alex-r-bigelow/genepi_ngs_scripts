#!/usr/bin/env python
import argparse, sys
from genome_utils import kgpInterface, countingDict

def run(args):
    kgp = kgpInterface(args.data,sys.path[0] + "/KGP_populations.txt")
    outfile = open(args.outfile,'w')
    freqOnly = args.frequencies_only.lower().startswith('t')
    
    wroteHeader = False
    for line,people in kgp.iterate():
        if not wroteHeader:
            outfile.write('CHROM\tPOS\tID')
            if not freqOnly:
                for p in people:
                    outfile.write('\t%s_1\t%s_2' % (p,p))
            outfile.write('\n')
            wroteHeader = True
        line.extractChrAndPos()
        line.extractAlleles()
        line.extractGenotypes()
        outfile.write([line.chromosome,str(line.position),line.name].join('\t'))
        if freqOnly:
            counts = countingDict()
            total = 0.0
            for p in people:
                if line.genotypes[p][0] != None:
                    counts[line.genotypes[p][0]] += 1
                    total += 1.0
                if line.genotypes[p][1] != None:
                    counts[line.genotypes[p][1]] += 1
                    total += 1.0
            for i,c in counts.iteritems():
                outfile.write('\t%s:\t%f' % (line.alleles[i],c/total))
        else:
            for p in people:
                a1 = line.genotypes[p][0]
                if a1 == None:
                    a1 = '.'
                else:
                    a1 = line.alleles[a1]
                a2 = line.genotypes[p][1]
                if a2 == None:
                    a2 = '.'
                else:
                    a2 = line.alleles[a2]
                
                outfile.write('\t%s\t%s' % (a1,a2))
        outfile.write('\n')
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates a .csv file containing all the genotypes for a population in the 1000 Genomes data set.')
    parser.add_argument('--out', type=str, dest="outfile", required=True,
                        help='output .csv file')
    parser.add_argument('--data', type=str, dest="data", required=True,
                        help='Path to directory containing 1000 Genomes .vcf.gz files 1-22,X,Y')
    parser.add_argument('--population', type=str, dest="pop", default="kALL",
                        help='One of the column headers in KGP_populations.txt. Defaults to kALL.')
    parser.add_argument('--frequencies_only', type=str, dest="frequencies_only", nargs="?", const="True", default="False",
                        help='Instead of dumping all the genotypes, just give a frequency for each allele.')
        
    args = parser.parse_args()
    run(args)