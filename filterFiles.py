#!/usr/bin/env python
import argparse, os, sys
from index_kgp import vcfLine, kgpInterface, countingDict

def tick():
    print ".",

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
    parser.add_argument('--expression', type=str, dest="expression",
                        help="File containing a Python-syntax expression to evaluate, e.g.:\n\n"+
                        "'%%s' == 'match this string' and %%s %% %%s == 0\n\n"+
                        "would evaluate to true (and the line would be included)\n"+
                        "if the first column in --columns matches 'match this string'\n"+
                        "exactly and the second column in --columns is an integer\n"+
                        "that can be evenly divided by the third column in --columns.\n"+
                        "Any valid python eval() code is permitted.")
    parser.add_argument('--columns', type=str, dest="columns",
                        help="File containing INFO field IDs or 'CHROM', 'POS', 'ID', 'QUAL' or 'FILTER' to\n"+
                        "use in --expression. Note that FILTER will always be a list\n"+
                        "(even if it just has one element), and INFO fields will be a\n"+
                        "list if a comma is present in the value (this could change from\n"+
                        "row to row). QUAL will be converted to a float, but all other\n"+
                        "values will be treated as strings. Any missing values will\n"+
                        "yield a string of a single period \".\" You'll need to consider\n"+
                        "conversions/error checking in your expressions.")
    
    args = parser.parse_args()
    
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    failfile = None
    if args.failfile != "":
        failfile = open(args.failfile,'w')
    errfile = None
    if args.errfile != "":
        errfile = open(args.errfile,'w')
    
    tempfile = open(args.expression,'r')
    expression = tempfile.readline()
    tempfile.close()
    
    columns = []
    tempfile = open(args.columns, 'r')
    for line in tempfile:
        line = line.strip()
        columns.append(line)
    tempfile.close()
    
    for line in infile:
        line = line.strip()
        if len(line) <= 1:
            continue
        elif line.startswith("##"):
            outfile.write()
        
        exp = expression % tuple([temp[i] for i in columns])
        try:
            result = eval(exp)
            if temp[0] == 'chr4' and result:
                print line
            if result == True:
                outfile.write(line)
            elif result == False:
                if failfile != None:
                    failfile.write(line)
            else:
                if errfile != None:
                    errfile.write(line)
        except:
            if errfile != None:
                errfile.write(line)
    
    infile.close()
    outfile.close()
    if errfile != None:
        errfile.close()