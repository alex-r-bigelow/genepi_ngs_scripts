#!/usr/bin/env python
import argparse, csv, gzip, os
from genome_utils import vcfLine, genomeException, standardizeChromosome, chromosomeOrder

class csvLine:
    def __init__(self, columns, chrColumn, posColumn, idColumn=None):
        self.columns = columns
        self.chrom = standardizeChromosome(columns[chrColumn])
        self.pos = int(columns[posColumn])
        if idColumn == None:
            self.id = "%s_%i" % (self.chrom,self.pos)
        else:
            self.id = "."

def sniffCsv(path):
    # first peek at the first 20 lines or so
    bloodhound = csv.Sniffer()
    snifflines = 20
    text = []
    if path.endswith('.gz'):
        csvfile = gzip.open(path)
    else:
        csvfile = open(path,'rb')
    for x in xrange(snifflines):
        text.append(csvfile.readline())
        if text[x] == '':
            if x == 0:
                raise Exception('Error parsing header from .csv file!')
            break
    
    scent = bloodhound.sniff("".join(text))
    
    headers = text[0].strip().split(scent.delimiter)
    if not 'CHROM' in headers or not 'POS' in headers:
        raise genomeException('Missing CHROM or POS header in .csv file!')
    chromColumn = headers.index('CHROM')
    posColumn = headers.index('POS')
    if 'ID' in headers:
        idColumn = headers.index('ID')
    else:
        idColumn = None
    csvfile.close()
    
    return (scent.delimiter,headers,chromColumn,posColumn,idColumn)

def run(args):
    delimiter,headers,chromColumn,posColumn,idColumn = sniffCsv(args.csvfile)
    if args.exact == None and args.nearest == None and args.interpolate == None:
        temp = headers
        temp.remove('CHROM')
        temp.remove('POS')
        args.exact = temp
    
    exact = {}
    nearest = {}
    interpolate = {}
    
    vcffile = open(args.infile,'r')
    csvfile = open(args.csvfile,'r')
    outfile = open(args.outfile,'w')
    
    csvbasename = os.path.split(args.csvfile)[1]
    
    takenTags = set()
    vcfHeaderLine = ""
    
    while True:
        line = vcffile.readline()
        if len(line) <= 1:
            continue
        elif line.startswith("##"):
            if line.startswith("##INFO"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                takenTags.add(newTag)
            outfile.write(line)
        elif line.startswith("#"):
            vcfHeaderLine = line
            break
        else:
            raise Exception("Missing a header line or something else is wrong...")
    
    if args.exact != None:
        for x in args.exact:
            if not x in headers:
                raise Exception('Column header "%s" doesn\'t exist in %s' % (x,csvbasename))
            temp = x
            dupNumber = 2
            while x in takenTags:
                x = temp + str(dupNumber)
                dupNumber += 1
            takenTags.add(x)
            outfile.write('##INFO=<ID=%s,Number=.,Type=String,Description=\"addCSVtoVCF.py: Column %s from %s in --exact match mode\">\n' % (x,temp,csvbasename))
            exact[x] = headers.index(x)
    
    if args.nearest != None:
        for x in args.nearest:
            if not x in headers:
                raise Exception('Column header "%s" doesn\'t exist in %s' % (x,args.csvfile))
            temp = x
            dupNumber = 2
            while x in takenTags:
                x = temp + str(dupNumber)
                dupNumber += 1
            takenTags.add(x)
            outfile.write('##INFO=<ID=%s,Number=.,Type=String,Description=\"addCSVtoVCF.py: Column %s from %s in --nearest match mode\">\n' % (x,temp,csvbasename))
            nearest[x] = headers.index(x)
    
    if args.interpolate != None:
        for x in args.interpolate:
            if not x in headers:
                raise Exception('Column header "%s" doesn\'t exist in %s' % (x,args.csvfile))
            temp = x
            dupNumber = 2
            while x in takenTags:
                x = temp + str(dupNumber)
                dupNumber += 1
            takenTags.add(x)
            outfile.write('##INFO=<ID=%s,Number=.,Type=String,Description=\"addCSVtoVCF.py: Column %s from %s in --interpolate match mode\">\n' % (x,temp,csvbasename))
            interpolate[x] = headers.index(x)
    
    outfile.write(vcfHeaderLine)
    
    # grab our first lines  
    vLine = vcffile.readline()
    vLine = vcfLine(vLine.strip().split('\t'))
    vLine.extractChrAndPos()
    vLine.extractInfo()
    
    lastCline = None
    csvfile.readline()  # skip the header
    cLine = csvLine(csvfile.readline().strip().split(delimiter),chromColumn,posColumn,idColumn)
    
    while True:
        speedAhead = cLine != None  # a flag that lets us just spit out .vcf lines because we know that either the .csv file has finished or there are no new .csv lines on the same chromosome
        # ... Are we even on the same chromosome?
        while speedAhead and cLine.chrom != vLine.chromosome:
            # the .csv file is ahead by a whole chromosome at least... keep lastCline intact and just spew out .vcf lines until it catches up
            if chromosomeOrder.index(cLine.chrom) > chromosomeOrder.index(vLine.chromosome):
                speedAhead = False
                break
            else:
                # okay, the .csv file is behind the .vcf file by at least a chromsome... speed ahead until we catch up or run out of .csv data
                lastCline = cLine
                cLine = csvfile.readline()
                if not cLine:
                    # shoot... we're out of .csv data. We already know that lastCline wasn't on the same chromosome as the current vLine, so make it None as well
                    cLine = None
                    lastCline = None
                    speedAhead = False
                    break
                else:
                    cLine = csvLine(cLine.strip().split(delimiter),chromColumn,posColumn,idColumn)
        # Okay, now we're on the same chromosome... zip ahead until cLine and lastCline are straddling vLine
        while speedAhead and cLine.pos < vLine.position:
            lastCline = cLine
            cLine = csvfile.readline()
            if not cLine:
                # shoot... out of .csv data. We know lastCline is still on the same chromosome, so preserve that, but make cLine None so we know nothing is left
                cLine = None
                break
            else:
                cLine = csvLine(cLine.strip().split(delimiter),chromColumn,posColumn,idColumn)
        
        # Whew! We're finally straddling the vLine...
        
        # Check the super-special case first (exact match)
        if cLine != None and cLine.pos == vLine.position:
            for x,i in exact.iteritems():
                vLine.info[x] = cLine.columns[i]
            for x,i in nearest.iteritems():
                vLine.info[x] = cLine.columns[i]
            for x,i in interpolate.iteritems():
                vLine.info[x] = cLine.columns[i]
        elif cLine != None: # cLine.pos will be > vLine.position
            if not args.omit_mismatches:
                for x,i in exact.iteritems():
                    vLine.info[x] = "."
            if lastCline == None:
                for x,i in nearest.iteritems():
                    vLine.info[x] = cLine.columns[i]
                for x,i in interpolate.iteritems():
                    vLine.info[x] = cLine.columns[i]
            else:
                closestLine = lastCline if vLine.position - lastCline.pos <= cLine.pos - vLine.position else cLine
                for x,i in nearest.iteritems():
                    vLine.info[x] = closestLine.columns[i]
                for x,i in interpolate.iteritems():
                    try:
                        lastVal = float(lastCline.columns[i])
                        nextVal = float(cLine.columns[i])
                        vLine.info[x] = str(lastVal + (nextVal - lastVal)*(vLine.position - lastCline.pos)/(cLine.pos - lastCline.pos))
                    except ValueError:
                        vLine.info[x] = closestLine.columns[i]
        else: # cLine == None
            if lastCline == None:
                if not args.omit_mismatches:
                    for x,i in exact.iteritems():
                        vLine.info[x] = "."
                    for x,i in nearest.iteritems():
                        vLine.info[x] = "."
                    for x,i in interpolate.iteritems():
                        vLine.info[x] = "."
            else:
                if not args.omit_mismatches:
                    for x,i in exact.iteritems():
                        vLine.info[x] = "."
                for x,i in nearest.iteritems():
                    vLine.info[x] = lastCline.columns[i]
                for x,i in interpolate.iteritems():
                    vLine.info[x] = lastCline.columns[i]
        # Okay, we've copied over everything; write the line
        outfile.write(str(vLine))
        # Grab the next one
        vLine = vcffile.readline()
        if not vLine:
            break   # No more variants - we're done!
        vLine = vcfLine(vLine.strip().split('\t'))
        vLine.extractChrAndPos()
        vLine.extractInfo()
    
    csvfile.close()
    vcffile.close()
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Attempts to incorporate .csv data as INFO fields in a .vcf file. Unless --exact, --nearest, or --interpolate are specified, '+
                                     'all columns will be incorporated in --exact mode. ALL FILES SHOULD HAVE PREVIOUSLY BEEN SORTED WITH sort.py')
    parser.add_argument('--in', type=str, dest="infile", required=True,
                        help='Input .vcf file.')
    parser.add_argument('--csv', type=str, dest="csvfile", required=True,
                        help='Input .csv file. There are some important things here that you, the .csv creator, should have '+
                                    'taken care of before trying to add your own stats. Python is pretty flexible '+
                                    'with figuring out delimiters on its own, but '+
                                    'we expect at least a few things: you need a single header line, with some '+
                                    'kind of label for each column. These labels will be used as INFO IDs as '+
                                    'closely as possible, but if the ID already exists, we\'ll start tacking on '+
                                    'numbers to the labels. You should also have at least a \"CHROM\" and a '+
                                    '\"POS\" column (case sensitive). You can optionally have an \"ID\" column '+
                                    'that will override existing rs numbers in the .vcf file. As we can\'t anticipate '+
                                    'how many values each field could have, by default the INFO Number= parameter will be '+
                                    '\".\" and the type will be \"String\". If this matters to you, it should be simple to change this yourself in the .vcf '+
                                    'header. Remember that multiple values in a single field need to be comma-delimited, so '+
                                    'you\'ll probably want to use something else like tabs for the rest of your file. '+
                                    'Your .csv file and .vcf file should both be sorted in the same chromosome order.')
    parser.add_argument('--out', type=str, dest="outfile", required=True,
                        help='Output .vcf file')
    parser.add_argument('--omit_mismatches', type=bool, dest="omit_mismatches", default=False,
                        help="Instead of marking variants that have no matches as missing \".\", omit the field entirely.")
    parser.add_argument('--exact', type=str, dest="exact", nargs="+",
                        help='Performs exact position matching for each column listed as an argument. A value will only be applied from the .csv file if a variant in the ' +
                            '.vcf file matches CHROM and POS exactly. Variants with no matching rows in the .csv file will be marked as missing that value ("."). Variants with ' +
                            'no ID in the .vcf file that have an exact match with an ID in the .csv file will be updated.')
    parser.add_argument('--nearest', type=str, dest="nearest", nargs="+",
                        help='Performs nearest neighbor matching for each column listed as an argument. Every variant in the .vcf file will be assigned the value of the nearest ' +
                            'matching .csv row. Variants on chromosomes that have no corresponding .csv lines will be marked as missing that value (".")')
    parser.add_argument('--interpolate', type=str, dest="interpolate", nargs="+",
                        help='Performs linear interpolation for each column listed as an argument. Every variant will be assigned an interpolated value from the nearest .csv lines ' +
                            'in either direction. Variants beyond .csv endpoints or .csv values that can\'t be converted to numbers will revert to --nearest matching behavior.')
    
    args = parser.parse_args()
    run(args)