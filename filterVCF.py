#!/usr/bin/env python
import argparse
from genome_utils import vcfLine, bedLine

def run(args):
    
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    failfile = None
    if args.failfile != "":
        failfile = open(args.failfile,'w')
    errfile = None
    if args.errfile != "":
        errfile = open(args.errfile,'w')
    
    if args.expression != "":
        tempfile = open(args.expression,'r')
        expression = tempfile.readline()
        tempfile.close()
    else:
        expression = "True"
    
    columns = []
    if args.columns != "":
        tempfile = open(args.columns, 'r')
        for line in tempfile:
            line = line.strip()
            columns.append(line)
        tempfile.close()
    
    bedRegions = None
    if args.bed != "":
        bedRegions = []
        tempfile = open(args.bed, 'r')
        for line in tempfile:
            bedRegions.append(bedLine(line.split()))
        tempfile.close()
    
    for line in infile:
        if len(line) <= 1:
            continue
        elif line.startswith("#"):
            outfile.write(line)
            if failfile != None:
                failfile.write(line)
            if errfile != None:
                errfile.write(line)
            continue
        else:
            line = vcfLine(line.strip().split('\t'))
            expArgs = []
            for c in columns:
                if c == "CHROM":
                    line.extractChrAndPos()
                    expArgs.append(line.chromosome)
                elif c == "POS":
                    line.extractChrAndPos()
                    expArgs.append(line.position)
                elif c == "ID":
                    line.extractChrAndPos()
                    expArgs.append(line.name)
                elif c == "QUAL":
                    line.extractQual()
                    expArgs.append(line.qual)
                elif c == "FILTER":
                    line.extractFilters()
                    expArgs.append(line.filters)
                else:
                    line.extractInfo()
                    expArgs.append(line.info.get(c,"."))
            
            # first see if it fails the .bed regions
            if bedRegions != None:
                passedBed = False
                line.extractChrAndPos()
                for bed in bedRegions:
                    if bed.contains(line.position):
                        passedBed = True
                        break
                if not passedBed:
                    if failfile != None:
                        failfile.write(str(line))
                        continue
            
            exp = expression % tuple(expArgs)
            try:
                result = eval(exp)
                if result == True:
                    outfile.write(str(line))
                elif result == False:
                    if failfile != None:
                        failfile.write(str(line))
                else:
                    if errfile != None:
                        errfile.write(str(line))
            except:
                if errfile != None:
                    errfile.write(str(line))
    
    infile.close()
    outfile.close()
    if errfile != None:
        errfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts a subset of a .vcf file using a custom expression.\n'+
                                                 'Note that the freedom of expression(a pun) that this program\n'+
                                                 'allows creates a SERIOUS SECURITY RISK for anything public.\n\n'+
                                                 'THIS PROGRAM SHOULD ONLY BE USED IN A DESKTOP ENVIRONMENT!', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in', type=str, dest="infile",
                        help='File to filter.')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='File to write lines that pass --expression')
    parser.add_argument('--fail', type=str, dest="failfile", nargs="?", default="",
                        help='File to write lines that fail --expression')
    parser.add_argument('--err', type=str, dest="errfile", nargs="?", default="",
                        help='File to write lines on which --expression generates an\nerror or produces a non-boolean result')
    parser.add_argument('--expression', type=str, dest="expression", nargs="?", default="",
                        help="File containing a Python-syntax expression to evaluate, e.g.:\n\n"+
                        "'%%s' == 'match this string' and %%s %% %%s == 0\n\n"+
                        "would evaluate to true (and the line would be included)\n"+
                        "if the first column in --columns matches 'match this string'\n"+
                        "exactly and the second column in --columns is an integer\n"+
                        "that can be evenly divided by the third column in --columns.\n"+
                        "Any valid python eval() code is permitted.")
    parser.add_argument('--columns', type=str, dest="columns", nargs="?", default="",
                        help="File containing INFO field IDs or 'CHROM', 'POS', 'ID', 'QUAL' or 'FILTER' to\n"+
                        "use in --expression. CHROM will always be converted to a string beginning\n"+
                        "with \"chr\", regardless of the .vcf format, POS will always be an integer,\n"+
                        "FILTER will always be a list (even if it just has one element), and INFO \n"+
                        "fields will be a list if a comma is present in the value, otherwise it will\n"+
                        "be a string (this could change from row to row). QUAL will be converted to a\n"+
                        "float. Any missing values will yield a string of a single period \".\"\n"+
                        "You'll need to consider conversions/error checking in your expressions.")
    parser.add_argument('--bed', type=str, dest="bed", nargs="?", default="",
                        help="In addition to --expression filters, this allows you to only include variants\n"+
                        "that lie within the regions specified in a .bed file.")
    
    args = parser.parse_args()
    run(args)