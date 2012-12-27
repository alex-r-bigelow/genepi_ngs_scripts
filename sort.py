#!/usr/bin/env python
import argparse, gzip, os, tempfile
from durus.file_storage import FileStorage
from durus.connection import Connection
from genome_utils import chromosomeOrder, standardizeChromosome, vcfLine
from blistPersistent.sorteddictPersistent import sorteddict
from addCSVtoVCF import sniffCsv

def sortVcf(inpath,outpath,compress,tickFunction=None,numTicks=1000):
    tickInterval = 2*os.path.getsize(inpath)/numTicks
    nextTick = 0
    
    if compress:
        infile = gzip.open(inpath,'rb')
        outfile = gzip.open(outpath,'wb')
    else:
        infile = open(inpath,'rb')
        outfile = open(outpath,'wb')
    
    tempFile = tempfile.NamedTemporaryFile()
    tempPath = tempFile.name
    tempFile.close()
    
    dataConnection = Connection(FileStorage(tempPath))
    data = dataConnection.get_root()
    extraChromosomes = []
    
    contigs = {None:[]}
    
    for line in infile:
        if line.startswith('##'):
            if 'contig' in line:
                contigID = line[line.find('ID=')+3:]
                before = line[:line.find('ID=')+3]
                after = contigID[contigID.find(','):]
                contigID = contigID[:contigID.find(',')]
                
                contigID = standardizeChromosome(contigID)
                
                if contigID in chromosomeOrder:
                    contigs[contigID] = before + contigID + after
                else:
                    contigs[None].append((contigID,before + contigID + after))
            else:
                outfile.write(line)
        elif line.startswith('#'):
            for c in chromosomeOrder:
                if contigs.has_key(c):
                    outfile.write(contigs[c])
            for c in sorted(contigs[None]):
                outfile.write(c[1])
            outfile.write(line)
        else:
            line = vcfLine(line.strip().split('\t'))
            line.extractChrAndPos()
            if not data.has_key(line.chromosome):
                data[line.chromosome] = sorteddict()
            data[line.chromosome][line.position] = line
            if infile.tell() > nextTick:
                dataConnection.commit()
                nextTick += tickInterval
                if tickFunction != None:
                    tickFunction()
    infile.close()
    dataConnection.commit()
    
    tickInterval = 0
    currentTick = 0
    nextTick = 0
    for lines in data.itervalues():
        tickInterval += len(lines)
    tickInterval = 2*tickInterval/numTicks
    
    for chrom in chromosomeOrder:
        if data.has_key(chrom):
            for line in data[chrom]:
                outfile.write(str(line))
                currentTick += 1
                if currentTick > nextTick:
                    nextTick += tickInterval
                    if tickFunction != None:
                        tickFunction()
    for chrom in extraChromosomes:
        if data.has_key(chrom):
            for line in data[chrom]:
                outfile.write(str(line))
                currentTick += 1
                if currentTick > nextTick:
                    nextTick += tickInterval
                    if tickFunction != None:
                        tickFunction()
    outfile.close()

def sortCsv(inpath,outpath,compress,tickFunction=None,numTicks=1000):
    tickInterval = 2*os.path.getsize(inpath)/numTicks
    nextTick = 0
    
    #
    
    if compress:
        infile = gzip.open(inpath,'rb')
        outfile = gzip.open(outpath,'wb')
    else:
        infile = open(inpath,'rb')
        outfile = open(outpath,'wb')
    
    tempFile = tempfile.NamedTemporaryFile()
    tempPath = tempFile.name
    tempFile.close()
    
    dataConnection = Connection(FileStorage(tempPath))
    data = dataConnection.get_root()
    extraChromosomes = []
    
    

def sortBed(inpath,outpath,compress,tickFunction=None,numTicks=1000):
    pass

def run(args):
    temp = os.path.splitext(args.infile)
    f = temp[1].lower()
    compress = False
    if f == ".gz":
        compress = True
        f = os.path.splitext(temp[0])[1].lower()
    if f == ".vcf":
        sortVcf(args.infile,args.outfile,compress)
    elif f == ".csv":
        sortCsv(args.infile,args.outfile,compress)
    elif f == ".bed":
        sortBed(args.infile,args.outfile,compress)
    else:
        raise Exception("Unknown format: %s" % f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sorts a .vcf, .csv, or .bed file first by chromosome: 1-22,X,Y,M, then by position (other chromosomes such as chrUn will be sorted alphabetically and placed last). If .csv, it should have \"CHROM\" and \"POS\" columns')
    parser.add_argument('--in', type=str, dest="infile",
                        help='Path to file (format automatically determined from extension). .gz compressed files are allowed, but your output will also be compressed')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='Path to file (output should be the same format as input). .gz compressed files are allowed, but your output will also be compressed')
    
    args = parser.parse_args()
    run(args)
    