#!/usr/bin/env python
import argparse, os
from index_kgp import vcfLine, kgpInterface, countingDict

def tick():
    print ".",

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
    print "Counting Alleles...",
    infile = open(args.infile,'r')
    for line in infile:
        columns = line.strip().split('\t')
        if len(columns) <= 1 or line.startswith("#"):
            continue
        line = vcfLine(columns)
        line.extractAlleles()
        numAlleles = max(numAlleles,len(line.alleles))
    infile.close()
    print numAlleles
    
    print "Opening KGP connection..."
    kgp = kgpInterface(args.data)
    
    print "Calculating..."
    outfile = open(args.outfile, 'w')
    
    outfile.write('Chromosome\tPosition')
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i'%i)
    for i in xrange(numAlleles):
        outfile.write('\t# Sharing Allele %i'%i)
    outfile.write('\t# Samples w/data')
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i Frequency'%i)
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i KGP Frequency'%i)
    
    total = 0
    noData = 0
    passedKGP = 0
    sharingCounts = countingDict()
    sharingFractions = {}
        
    for line,kgpLine in kgp.iterateVcf(args.infile,tickFunction=tick,numTicks=1000):
        outfile.write('\n')
        
        # write chromosome, position
        outfile.write('%s\t%i'%(line.chromosome,line.position))
        
        line.extractAlleles()
        
        for a in line.alleles:
            outfile.write('\t%s'% a)
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t-')
        
        # write sharing per allele
        sharing = line.getSharing()
        for a in line.alleles:
            outfile.write('\t%i'%sharing[a])
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t-')
        
        # write sharing denominator
        numberWithData = line.numberOfSamplesWithData()
        outfile.write('\t%i'%numberWithData)
        
        # write allele frequencies in file
        frequencies = line.getAlleleFrequencies()
        for a in line.alleles:
            outfile.write('\t%f'%frequencies[a])
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t-')
        
        # write kgp allele frequencies
        kgpFrequencies = None
        if kgpLine == None:
            for a in line.alleles:
                outfile.write('\tND')
            for x in xrange(numAlleles-len(line.alleles)):
                outfile.write('\t-')
        else:
            kgpLine.extractGenotypes()
            kgpFrequencies = kgpLine.getAlleleFrequencies()
            for a in line.alleles:
                if kgpFrequencies.has_key(a):
                    outfile.write('\t%f'%kgpFrequencies[a])
                else:
                    outfile.write('\tND')
            for x in xrange(numAlleles-len(line.alleles)):
                outfile.write('\t-')
        
        # count stuff for final output
        total += 1
        if kgpLine == None or numberWithData == 0:
            noData += 1
        else:
            lowestFreq = 1.0
            minorAllele = None
            for a in line.alleles:
                if kgpFrequencies.has_key(a):
                    if kgpFrequencies[a] < lowestFreq:
                        lowestFreq = kgpFrequencies[a]
                        minorAllele = a
            if lowestFreq <= 0.01:
                passedKGP += 1
                
                proportion = sharing[minorAllele]/float(numberWithData)
                sharingCounts[proportion] += 1
                sharingFractions[proportion] = (sharing[minorAllele],numberWithData)
    
    infile.close()
    outfile.close()
    print ""
    print "Done"
    print ""
    print "Total number of variants: %i" % total
    print "# of variants with at least some data from both sources: %i" % total-noData
    print "# left if variants with no allele <= 0.01 are removed: %i" % passedKGP
    temptotal = 0
    for p,c in sorted(sharingCounts.iteritems()):
        print "# left if variants sharing < %i/%i of the minor KGP allele are removed: %i" % (sharingFractions[p][0],sharingFractions[p][1],passedKGP-temptotal)
        temptotal += c