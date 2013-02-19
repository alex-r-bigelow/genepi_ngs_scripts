#!/usr/bin/env python
import argparse, os
from genome_utils import standardizeChromosome, vcfLine, bedLine, chromosomeRank
from addCSVtoVCF import sniffCsv
from recipe576755 import batch_sort

count = 0
def tick():
    global count
    print "Sorting: {0}%\r".format(count),
    count += 1
    return True

class vcfKey:
    FIRSTLINE = 0
    CONTIG = 1
    OTHER_META = 2
    HEADER = 3
    REGULAR = 4
    EMPTY = 5
    
    def __init__(self, line):
        if len(line.strip()) == 0:
            self.type = vcfKey.EMPTY
        elif line.startswith('##'):
            line = line.lower()
            if line.startswith('##fileformat'):
                self.type = vcfKey.FIRSTLINE
            elif line.startswith('##contig'):
                self.type = vcfKey.CONTIG
                contigID = line[line.find('ID=')+3:]
                contigID = contigID[:contigID.find(',')]
                self.chromosome = standardizeChromosome(contigID)
            else:
                self.type = vcfKey.OTHER_META
                self.line = line
        elif line.startswith('#'):
            self.type = vcfKey.HEADER
        else:
            self.type = vcfKey.REGULAR
            temp = vcfLine(line.strip().split('\t'))
            temp.extractChrAndPos()
            self.chromosome = temp.chromosome
            self.position = temp.position
    def __cmp__(self, other):
        if self.type != other.type:
            return self.type-other.type
        else:
            if self.type == vcfKey.REGULAR:
                if self.chromosome == other.chromosome:
                    return self.position-other.position
                else:
                    myRank = chromosomeRank.get(self.chromosome,len(chromosomeRank))
                    otherRank = chromosomeRank.get(other.chromosome,len(chromosomeRank))
                    result = myRank-otherRank
                    if result == 0:
                        # we already know they're not the same chromosome; our behavior is to sort unknown chromosomes alphabetically
                        return cmp(self.chromosome,other.chromosome)
                    else:
                        return result
            elif self.type == vcfKey.CONTIG:
                myRank = chromosomeRank.get(self.chromosome,len(chromosomeRank))
                otherRank = chromosomeRank.get(other.chromosome,len(chromosomeRank))
                result = myRank-otherRank
                if result == 0:
                    # we already know they're not the same chromosome; our behavior is to sort unknown chromosomes alphabetically
                    return cmp(self.chromosome,other.chromosome)
                else:
                    return result
            elif self.type == vcfKey.OTHER_META:
                # I could probably do a better job of this, but the .vcf spec doesn't require it so I don't care
                return cmp(self.line,other.line)
            elif self.type == vcfKey.EMPTY:
                return 0
            else:
                raise Exception("Duplicate ##fileformat or header lines!")

class csvKey:
    delimiter = None
    chromColumn = None
    posColumn = None
    
    def __init__(self, line):
        if "CHROM" in line and "POS" in line:
            self.isFirstLine = True
            self.line = line
        else:
            self.isFirstLine = False
            columns = line.strip().split(csvKey.delimiter)
            self.chromosome = standardizeChromosome(columns[csvKey.chromColumn])
            self.position = int(columns[csvKey.posColumn])
    
    def __cmp__(self, other):
        if self.isFirstLine:
            return -1
        elif other.isFirstLine:
            return 1
        elif self.chromosome == other.chromosome:
            return self.position-other.position
        else:
            myRank = chromosomeRank.get(self.chromosome,len(chromosomeRank))
            otherRank = chromosomeRank.get(other.chromosome,len(chromosomeRank))
            result = myRank-otherRank
            if result == 0:
                # we already know they're not the same chromosome; our behavior is to sort unknown chromosomes alphabetically
                return cmp(self.chromosome,other.chromosome)
            else:
                return result

class bedKey:
    def __init__(self, line):
        temp = bedLine(line)
        self.chromosome = temp.chromosome
        self.start = temp.start
    
    def __cmp__(self, other):
        if self.chromosome == other.chromosome:
            return self.start - other.start
        else:
            myRank = chromosomeRank.get(self.chromosome,len(chromosomeRank))
            otherRank = chromosomeRank.get(other.chromosome,len(chromosomeRank))
            result = myRank-otherRank
            if result == 0:
                # we already know they're not the same chromosome; our behavior is to sort unknown chromosomes alphabetically
                return cmp(self.chromosome,other.chromosome)
            else:
                return result

def sortVcf(inpath, outpath, tickFunction=tick, numTicks=100):
    batch_sort(inpath, outpath, key=vcfKey, tickFunction=tickFunction, numTicks=numTicks)

def sortCsv(inpath, outpath, tickFunction=tick, numTicks=100):
    csvKey.delimiter,headers,csvKey.chromColumn,csvKey.posColumn,idColumn = sniffCsv(inpath)
    batch_sort(inpath, outpath, key=csvKey, tickFunction=tickFunction, numTicks=numTicks)

def sortBed(inpath, outpath, tickFunction=tick, numTicks=100):
    batch_sort(inpath, outpath, key=bedKey, tickFunction=tickFunction, numTicks=numTicks)

def run(args, tickFunction=tick, numTicks=100):
    inpath = args.infile
    outpath = args.outfile
    
    temp = os.path.splitext(inpath)
    f = temp[1].lower()
    if f == ".vcf":
        sortVcf(inpath,outpath,tickFunction,numTicks)
    elif f == ".csv":
        sortCsv(inpath,outpath,tickFunction,numTicks)
    elif f == ".bed":
        sortBed(inpath,outpath,tickFunction,numTicks)
    else:
        raise Exception("Unknown format: %s" % f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sorts a .vcf, .csv, or .bed file first by chromosome: 1-22,X,Y,M, then by position (other chromosomes such as chrUn will be sorted alphabetically and placed last). If .csv, it should have \"CHROM\" and \"POS\" columns')
    parser.add_argument('--in', type=str, dest="infile", required=True,
                        help='Path to file (format automatically determined from extension)')
    parser.add_argument('--out', type=str, dest="outfile", required=True,
                        help='Path to file (output should be the same format as input)')
    
    args = parser.parse_args()
    run(args)
    