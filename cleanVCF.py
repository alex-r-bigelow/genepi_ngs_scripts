#!/usr/bin/env python
import argparse, gzip, os
from genome_utils import vcfLine, infoDetails, MAX_INFO_STRINGS

def extractInfoFields(path,max_strings=MAX_INFO_STRINGS,tickFunction=None,numTicks=1000):
    infoFields = {}
    
    tickInterval = os.path.getsize(path)/numTicks
    nextTick = 0
    
    if path.endswith('gz'):
        infile = gzip.open(path,'rb')
    else:
        infile = open(path,'rb')
    for line in infile:
        if tickFunction != None:
            if infile.tell() > nextTick:
                if not tickFunction():
                    infile.close()
                    return None
                nextTick += tickInterval
        line = line.strip()
        if len(line) <= 1:
            continue
        elif line.startswith("##INFO"):
            newTag = line[line.find("ID=")+3:]
            newTag = newTag[:newTag.find(',')]
            if infoFields.has_key(newTag):
                raise Exception("Duplicate INFO ID or use of reserved ID:\t%s" % newTag)
            infoFields[newTag] = infoDetails(newTag, max_strings, False)
        elif line.startswith("#"):
            continue
        else:
            line = vcfLine(line.split('\t'))
            line.extractInfo()
            for k,v in line.info.iteritems():
                if not infoFields.has_key(k):
                    raise Exception("Missing ##INFO pragma for: %s" % k)
                else:
                    infoFields[k].addArbitraryValue(v)
    infile.close()
    return infoFields

def run(args):
    print 'Counting values...'
    max_strings = args.max_strings
    ignoreStringCounts = False
    if max_strings <= 0:
        ignoreStringCounts = True
        max_strings = MAX_INFO_STRINGS
    infoFields = extractInfoFields(args.infile, max_strings)
    
    validFields = set()
    for k,f in infoFields.iteritems():
        if not ignoreStringCounts and f.maxedOut:
            continue
        if args.preserve_info != None and k not in args.preserve_info:
            continue
        if args.remove_info != None and k in args.remove_info:
            continue
        validFields.add(k)
    
    print 'Writing file...'
    outfile = open(args.outfile, 'w')
    infile = open(args.infile, 'r')
    for line in infile:
        line = line
        if len(line) <= 1:
            continue
        elif line.startswith("##INFO"):
            newTag = line[line.find("ID=")+3:]
            newTag = newTag[:newTag.find(',')]
            if not infoFields.has_key(newTag):
                raise Exception("Second pass lost info tag:\t%s" % newTag)
            if newTag in validFields:
                outfile.write(line)
        elif line.startswith("#"):
            outfile.write(line)
        else:
            line = vcfLine(line.strip().split('\t'))
            line.extractInfo()
            keys = line.info.keys()
            for k in keys:
                if k not in validFields:
                    del line.info[k]
            outfile.write(str(line))
    infile.close()
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Cleans a .vcf file for viewing in an early version of compreheNGSive (0.2.2 - automagically removes categorical info fields with an excessive number of possible values)')
    parser.add_argument('--in', type=str, dest="infile", required = True,
                        help='Path to .vcf file')
    parser.add_argument('--out', type=str, dest="outfile", required = True,
                        help='Path to .vcf file')
    parser.add_argument('--max_strings', type=int, dest="max_strings", nargs="?", const=MAX_INFO_STRINGS, default=MAX_INFO_STRINGS,
                        help='Maximum number of strings a categorical INFO field can have before it\'s removed. If zero or negative, no limit is enforced. Default is %i' % MAX_INFO_STRINGS)
    parser.add_argument('--remove_info', type=str, dest="remove_info", nargs="+",
                        help='Remove specific field(s) from the .vcf file')
    parser.add_argument('--preserve_info', type=str, dest="preserve_info", nargs="+",
                        help='Only include the specified field(s) from the .vcf file (though they still can be excluded by --max_strings or --remove_info)')
    
    args = parser.parse_args()
    run(args)
    