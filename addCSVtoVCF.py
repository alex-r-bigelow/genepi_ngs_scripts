#!/usr/bin/env python
import argparse, os, csv
from index_kgp import vcfLine, standardizeChromosome

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Attempts to incorporate .csv data as INFO fields in a .vcf file. '+
                                    'There are some important things here that you, the .csv creator, should have '+
                                    'taken care of before trying to add your own stats. Python is pretty flexible '+
                                    'with figuring out delimiters on its own, but '+
                                    'we expect at least a few things: you need a single header line, with some '+
                                    'kind of label for each column. These labels will be used as INFO IDs as '+
                                    'closely as possible, but if the ID already exists, we\'ll start tacking on '+
                                    'numbers to the labels. You should also have at least a \"CHROM\" and a '+
                                    '\"POS\" column (case sensitive). Note that exact matches are expected; '+
                                    'if an indel position is one base pair off, that .csv data won\'t be added. You can optionally have an \"ID\" column '+
                                    'that will override existing rs numbers in the .vcf file. As we can\'t anticipate '+
                                    'how many values each field could have, by default the INFO Number= parameter will be '+
                                    '\".\" If this matters to you, it should be simple to change this yourself in the .vcf '+
                                    'header. Remember that multiple values in a single field need to be comma-delimited, so '+
                                    'you\'ll probably want to use something else like tabs for the rest of your file. '+
                                    'Your .csv file and .vcf file should both be sorted in the same chromosome order.')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--csv', type=str, dest="csvfile",
                        help='input .csv file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    
    args = parser.parse_args()
    
    # first peek at the first 20 lines or so
    bloodhound = csv.Sniffer()
    snifflines = 20
    text = []
    csvfile = open(args.csvfile,'r')
    for x in xrange(snifflines):
        text.append(csvfile.readline())
        if text[x] == '':
            raise Exception('Less than %i lines in .csv file' % snifflines)
    
    scent = bloodhound.sniff("".join(text))
    
    headers = text[0].strip().split(scent.delimiter)
    chromColumn = headers.index('CHROM')
    posColumn = headers.index('POS')
    if 'ID' in headers:
        idColumn = headers.index('ID')
    else:
        idColumn = None
    
    # now I have to do a pass through the .csv to pull out its chromosome order (we're assuming the .csv is smaller than the .vcf)
    chromOrder = []
    allChroms = set()
    for t in text[1:]:
        chrom = standardizeChromosome(t.strip().split(scent.delimiter)[chromColumn])
        if chrom not in allChroms:
            allChroms.add(chrom)
            chromOrder.append(chrom)
    
    for t in csvfile:
        chrom = standardizeChromosome(t.strip().split(scent.delimiter)[chromColumn])
        if chrom not in allChroms:
            allChroms.add(chrom)
            chromOrder.append(chrom)
    csvfile.close()
    
    # Okay, on to the main pass through both files
    
    takenTags = set()
    
    csvfile = open(args.csvfile,'r')
    vcffile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    
    cLine = csvfile.readline()  # skip the header
    cLine = csvfile.readline().strip().split(scent.delimiter)
    chrom = standardizeChromosome(cLine[chromColumn])
    pos = int(cLine[posColumn])
    eof = False
    
    for line in vcffile:
        if len(line) <= 1:
            continue
        elif line.startswith('##'):
            if line.startswith("##INFO"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                takenTags.add(newTag)
            outfile.write(line)
        elif line.startswith('#'):
            for h in headers:
                i = 2
                newTag = h
                while newTag in takenTags:
                    newTag = h + str(i)
                    i+=1
                outfile.write("##INFO=<ID=%s,Number=.,Type=Float,Description=\"User column added with addCSVtoVCF.py\">\n" % newTag)
            outfile.write(line)
        else:
            if eof: # we've exhausted the .csv file; just write what's left in the .vcf
                outfile.write(line)
                continue
            
            line = vcfLine(line.strip().split('\t'))
            line.extractChrAndPos()
            
            if line.chromosome != chrom:
                if line.chromosome not in allChroms:
                    # the .csv file is ahead; just output the .vcf line
                    outfile.write(str(line))
                    continue
                else:
                    # the .vcf file is ahead; output the .vcf line
                    outfile.write(str(line))
                    
                    # throw out all the .csv chromosomes before and including this one
                    index = 1
                    for c in chromOrder:
                        index += 1
                        allChroms.remove(c)
                    chromOrder = chromOrder[index:]
                    if len(chromOrder) == 0:
                        eof = True
                        continue
                    
                    # fast forward until we get to the new chromOrder[0] or eof
                    while chrom != chromOrder[0]:
                        text = csvfile.readline()
                        if text == '':
                            eof = True
                            break
                        cLine = text.strip().split(scent.delimiter)
                        chrom = standardizeChromosome(cLine[chromColumn])
                    if eof:
                        continue
                    pos = int(cLine[posColumn])
                    continue
            
            # okay, we're on the same chromosome; fast forward until we reach or pass the .vcf line
            wroteLine = False
            while pos < line.position:
                text = csvfile.readline()
                if text == '':
                    eof = True
                    break
                cLine = text.strip().split(scent.delimiter)
                chrom = standardizeChromosome(cLine[chromColumn])
                pos = int(cLine[posColumn])
                if chrom != line.chromosome:
                    # whoops, we just ran past the chromosome altogether; just output the .vcf line and move along
                    wroteLine = True
                    outfile.write(str(line))
                    break
            if wroteLine:
                continue
            
            # okay, now we're as close as we could possibly be; did we match?
            if pos == line.position:
                # A hit! A palpable hit!
                if idColumn != None:
                    line.id = cLine[idColumn]   # override the id if that's what we're doing
                
                line.extractInfo()
                for i,c in enumerate(cLine):
                    if i == chromColumn or i == posColumn or i == idColumn:
                        continue
                    if "," in c:
                        c = c.split(',')
                    line.info[headers[i]] = c
            
            # okay, if the .vcf line has a match it's been updated, so just write it and move along
            outfile.write(str(line))
            
    vcffile.close()
    csvfile.close()
    outfile.close()