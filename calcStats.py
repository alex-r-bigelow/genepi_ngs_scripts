#!/usr/bin/env python
import argparse
from index_kgp import vcfLine, kgpInterface

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a .csv file with some calculated statistics.')
    parser.add_argument('--data', type=str, dest="data",
                        help='Directory housing decompressed .vcf files 1-22,X. You should have created an index by running index_kgp before running this.')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .csv file')
    
    args = parser.parse_args()
    
    numAlleles = 0
    print "Counting Alleles..."
    infile = open(args.infile,'r')
    for line in infile:
        line = vcfLine(line)
        line.extractAlleles()
        numAlleles = max(numAlleles,len(line.alleles))
    infile.close()
    
    print "Opening KGP connection..."
    kgp = kgpInterface(args.data)
    
    print "Calculating..."
    infile = open(args.infile,'r')
    outfile = open(args.outfile, 'w')
    
    outfile.write('Chromosome\tPosition')
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i'%i)
    outfile.write('\t# Samples w/data')
    for i in xrange(numAlleles):
        outfile.write('\t# Sharing Allele %i'%i)
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i Frequency'%i)
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i KGP Frequency'%i)
    
    for line in infile:
        if len(line.strip()) <= 1 or line.startswith("#"):
            continue
        outfile.write('\n')
        
        line = vcfLine(line)
        line.extractChrAndPos()
        outfile.write('%s\t%i'%(line.chromosome,line.position))
        
        line.extractAlleles()
        outfile.write('\t%s'%('\t'.join(line.alleles)))
        
        outfile.write('\t%i'%line.numberOfSamplesWithData())
        
        line.extractGenotypes()
        sharing = line.getSharing(line.genotypes.iterkeys())
        for a in line.alleles:
            outfile.write('\t%i'%sharing[a])
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t')
        
        frequencies = line.getAlleleFrequencies(line.genotypes.iterkeys())
        for a in line.alleles:
            outfile.write('\t%f'%frequencies[a])
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t')
        
        kgpLines = kgp.getVcfLines(line.chromosome,line.position,line.position+1)
        if len(kgpLines) != 1:
            if len(kgpLines) == 0:
                for x in xrange(numAlleles):
                    outfile.write('\tND')
            else:
                raiseException('More than one KGP line came back:\n%s' % (str(kgpLines)))
        else:
            kgpLines[0].extractGenotypes()
            kgpFrequencies = kgpLines[0].getAlleleFrequencies(kgpLines[0].genotypes.iterkeys())
            for a in line.alleles:
                if kgpFrequencies.has_key(a):
                    outfile.write('\t%f'%kgpFrequencies[a])
                else:
                    outfile.write('\tNA')
            for x in xrange(numAlleles-len(line.alleles)):
                outfile.write('\t')
    infile.close()
    outfile.close()
    print "Done"