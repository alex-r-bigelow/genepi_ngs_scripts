#!/usr/bin/env python
import argparse, os, datetime
from index_kgp import vcfLine, kgpInterface, countingDict

def tick():
    print ".",

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a .cvf file from the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO fields of a .vcf file')
    parser.add_argument('--in', type=str, dest="infile",
                        help='Path to .vcf file')
    parser.add_argument('--out', type=str, dest="outfiles", nargs="+",
                        help='Path to .cvf file')
    
    args = parser.parse_args()
    
    numColumns = {"ALT":0,"FILTER":0}
    
    # TODO: get the numeric ranges, all valid categorical values
    
    print "Counting columns, ranges, values...",
    infile = open(args.infile,'r')
    for line in infile:
        line = line.strip()
        if len(columns) <= 1:
            continue
        elif line.startswith("#"):
            if line.startswith("##INFO"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                if numColumns.has_key(newTag):
                    raise Exception("Duplicate INFO ID or use of reserved ID:\t%s" % newTag)
                numColumns[newTag] = 0
            continue
        else:
            line = vcfLine(line.split('\t'))
            line.extractAlleles()
            line.extractInfo()
            line.extractFilters()
            numColumns['ALT'] = max(numColumns['ALT'],len(line.alleles))
            numColumns['FILTER'] = max(numColumns,len(line.filters))
            for k,v in numColumns.iteritems():
                if k == 'ALT':
                    numColumns[k] = max(v,len(line.alleles))
                elif k == 'FILTER':
                    numColumns[k] = max(v,len(line.filters))
                else:
                    if not numColumns.has_key(k):
                        raise Exception("Missing ##INFO pragma for ID: %s" % k)
                    elif isinstance(line.info[k],list):
                        numColumns[k] = max(v,len(line.info[k]))
                    elif v != None:
                        numColumns[k] = max(v,1)
    infile.close()
    
    print "Creating file..."
    outfile = open(args.outfile, 'w')
    
    outfile.write("##\t%s created from %s on %s\n" % (args.outfile,args.str(infile,datetime.datetime.now())))
    outfile.write("#\tChromosome\tCHR\n")
    outfile.write("#\tPosition\tPOS\n")
    outfile.write("#\tID\tID\n")
    for i in xrange(len(numColumns['ALT'])):
        outfile.write("#\tAllele %i\tIGNORE\n" % (i+1))
    outfile.write("#\tQUAL\tNUMERIC\t0.0\t%f\n" % maxQual)
    for k,v in sorted(numColumns.iteritems()):
        
    
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
    print "# of variants with at least some data from both sources: %i" % (total-noData)
    print "# left if variants with no allele <= 0.01 are removed: %i" % passedKGP
    temptotal = 0
    for p,c in sorted(sharingCounts.iteritems()):
        print "# left if variants sharing < %i/%i of the minor KGP allele are removed: %i" % (sharingFractions[p][0],sharingFractions[p][1],passedKGP-temptotal)
        temptotal += c