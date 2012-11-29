#!/usr/bin/env python
import argparse, os
from index_kgp import vcfLine, kgpInterface, countingDict

def tick():
    print ".",

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add some calculated statistics to a .vcf file\'s INFO column.')
    parser.add_argument('--data', type=str, dest="data",
                        help='Directory housing decompressed .vcf files 1-22,X. You should have created an index by running index_kgp before running this.')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    
    args = parser.parse_args()
    
    takenTags = set()
    headerline = ""
    
    print "Parsing Header...",
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    for line in infile:
        columns = line.strip().split('\t')
        if len(columns) <= 1:
            continue
        elif line.startswith("##"):
            if line.startswith("##INFO"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                takenTags.add(newTag)
            outfile.write(line)
            continue
        elif line.startswith("#"):
            headerline = line
            break
        else:
            raise Exception("Missing a header line or something else is wrong...")
    infile.close()
    
    kgpTag = "1kgAF"
    i = 2
    while kgpTag in takenTags:
        kgpTag = "1kgAF" + str(i)
        i += 1
    outfile.write("##INFO=<ID=%s,Number=.,Type=Float,Description=\"Frequency of each allele in the 1000 Genomes Project (calcStats.py)\">\n" % kgpTag)
    
    wDataTag = "wData"
    i = 2
    while wDataTag in takenTags:
        wDataTag = "wData" + str(i)
        i += 1
    outfile.write("##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Number of samples with calls (calcStats.py)\">\n" % wDataTag)
    
    sharingTag = "Sharing"
    i = 2
    while sharingTag in takenTags:
        sharingTag = "Sharing" + str(i)
        i += 1
    outfile.write("##INFO=<ID=%s,Number=.,Type=Integer,Description=\"Assuming all samples in this file have a common ancestor, the maximum possible number of samples that could share this allele (calcStats.py)\">\n" % sharingTag)
    
    outfile.write(headerline)
    
    print "Opening KGP connection..."
    kgp = kgpInterface(args.data)
    
    print "Calculating..."
    for line,kgpLine in kgp.iterateVcf(args.infile,tickFunction=tick,numTicks=1000):
        line.extractAlleles()
        line.extractInfo()
        
        # calculate sharing per allele
        line.info[sharingTag] = []
        sharing = line.getSharing()
        for a in line.alleles:
            line.info[sharingTag].append(str(sharing[a]))
        
        # write sharing denominator
        numberWithData = line.numberOfSamplesWithData()
        line.info[wDataTag] = str(numberWithData)
        
        # write kgp allele frequencies
        kgpFrequencies = None
        if kgpLine == None:
            line.info[kgpTag] = ['.' for a in line.alleles]
        else:
            kgpLine.extractGenotypes()
            kgpFrequencies = kgpLine.getAlleleFrequencies()
            line.info[kgpTag] = []
            for a in line.alleles:
                if kgpFrequencies.has_key(a):
                    line.info[kgpTag].append(str(kgpFrequencies[a]))
                else:
                    line.info[kgpTag].append('.')
        
        outfile.write(str(line)+'\n')
    
    infile.close()
    outfile.close()
    print ""
    print "Done"