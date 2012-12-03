#!/usr/bin/env python
import argparse
from genome_utils import vcfLine, bedLine

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Attempts to incorporate .bed data as INFO fields in a .vcf file; '+
                                    'if a variant resides in a .bed region, the score will be assigned to the variant '+
                                    'with a tag as close to the .bed region name (INFO IDs already in the .vcf file will '+
                                    'result in modified tags). As you could embed multiple values separated by commas, '+
                                    'the INFO Number= parameter will be \".\" If you care about this, you\'ll need to '+
                                    'change this yourself manually.')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--bed', type=str, dest="bedfile",
                        help='input .bed file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    
    args = parser.parse_args()
    
    scoreNames = set()
    bedRegions = []
    tempfile = open(args.bedfile, 'r')
    for line in tempfile:
        line = bedLine(line.split())
        bedRegions.append(line)
        scoreNames.add(line.name)
    tempfile.close()
    
    vcffile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    
    takenTags = set()
    
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
            for n in scoreNames:
                dupCount = 2
                newTag = n
                while newTag in takenTags:
                    newTag = n + str(dupCount)
                    dupCount += 1
                takenTags.add(newTag)
                outfile.write("##INFO=<ID=%s,Number=.,Type=Float,Description=\"User column added with addBEDtoVCF.py\">\n" % newTag)
            outfile.write(line)
        else:
            line = vcfLine(line.strip().split('\t'))
            line.extractChrAndPos()
            
            for b in bedRegions:
                if b.contains(line.position):
                    line.extractInfo()
                    line.info[b.name] = str(b.score)
            
            outfile.write(str(line))
            
    vcffile.close()
    outfile.close()