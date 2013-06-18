#!/usr/bin/env python
import argparse, sys
from genome_utils import standardizeChromosome, vcfLine, infoDetails

def run(args):
    separateInfoFields = args.separate_info_fields.strip().lower() == "true"
    numberAlleles = args.numbered_alleles.strip().lower() == "true"
    includeGenotypes = args.include_genotypes.strip().lower() == "true"
    includeGenotypeAttributes = args.include_genotype_attributes.strip().lower() == "true"
    if includeGenotypeAttributes:
        includeGenotypes = True
    
    ignoreFields = args.ignore_fields
    if ignoreFields == None:
        ignoreFields = set([])
    ignoreFields = set(ignoreFields)
    
    numAltAlleles = 1
    numFilters = 1
    infoOrder = []
    infoHeaders = {}
    peopleOrder = []
    formatOrder = []
    formatHeaders = {}
    
    infile = open(args.infile,'r')
    for line in infile:
        line = line.strip()
        if len(line) <= 1:
            continue
        elif line.startswith("#"):
            temp = line.lower()
            if temp.startswith("##info"):
                newTag = line[temp.find("id=")+3:]
                newTag = newTag[:newTag.find(',')]
                if not newTag in ignoreFields:
                    if not infoHeaders.has_key(newTag):
                        infoOrder.append(newTag)
                    infoHeaders[newTag] = 1
            elif temp.startswith("##format"):
                newTag = line[temp.find("id=")+3:]
                newTag = newTag[:newTag.find(',')]
                if not newTag in ignoreFields and not newTag == 'GT':
                    if not formatHeaders.has_key(newTag):
                        formatOrder.append(newTag)
                        formatHeaders[newTag] = 1
            elif temp.startswith("#chrom"):
                peopleOrder = line.split('\t')[9:]
        else:
            line = vcfLine(line.split('\t'))
            
            if separateInfoFields:
                line.extractAlleles()
                numAltAlleles = max(numAltAlleles,len(line.alleles)-1)  # don't include the REF allele
                
                line.extractFilters()
                numFilters = max(numFilters,len(line.filters))
            
            line.extractInfo()
            for k,v in line.info.iteritems():
                if not infoHeaders.has_key(k):
                    infoOrder.append(newTag)
                    infoHeaders[k] = 1
                if separateInfoFields and isinstance(v,list):
                    infoHeaders[k] = max(infoHeaders[k],len(v))
            
            if includeGenotypeAttributes:
                line.extractFormat()
                line.extractGenotypes()
                for i,p in enumerate(peopleOrder):
                    allele0,allele1,phased,attrs = line.genotypes[i]  # @UnusedVariable
                    for j,f in enumerate(line.format[1:]):
                        if len(attrs) > j:
                            formatHeaders[f] = max(formatHeaders[f],len(attrs[j].split(',')))
    infile.close()
    
    print "Creating file..."
    outfile = open(args.outfile, 'w')
    outfile.write('Chromosome\tPosition\tID\tReference_Allele')
    if separateInfoFields and numAltAlleles > 1:
        for x in xrange(numAltAlleles):
            outfile.write('\tAlternate_Allele_%i' % (x+1))
    else:
        outfile.write('\tAlternate_Allele')
    
    outfile.write('\tQual')
    
    if separateInfoFields and numFilters > 1:
        for x in xrange(numFilters):
            outfile.write('\tFilter_%i' % (x+1))
    else:
        outfile.write('\tFilter')
    
    for i in infoOrder:
        if separateInfoFields and infoHeaders[i] > 1:
            for x in xrange(infoHeaders[i]):
                outfile.write('\t%s_%i' % (i,x+1))
        else:
            outfile.write('\t%s' % i)
    
    if includeGenotypes:
        for p in peopleOrder:
            outfile.write('\t%s_Allele_1' % p)
            outfile.write('\t%s_Allele_2' % p)
            if includeGenotypeAttributes:
                outfile.write('\t%s_Phased' % p)
                for f in formatOrder:
                    if separateInfoFields and formatHeaders[f] > 1:
                        for x in xrange(formatHeaders[f]):
                            outfile.write('\t%s_%s_%i' % (p,f,x+1))
                    else:
                        outfile.write('\t%s_%s' % (p,f))
    outfile.write('\n')
    
    infile = open(args.infile,'r')
    for line in infile:
        line = line.strip()
        if len(line) <= 1 or line.startswith("#"):
            continue
        line = vcfLine(line.split('\t'))
        
        line.extractChrAndPos()
        outfile.write("%s\t%i\t%s" % (line.chromosome,line.position,line.id))
        
        line.extractAlleles()
        outfile.write('\t%s' % line.alleles[0])
        if separateInfoFields:
            for a in line.alleles[1:]:
                outfile.write('\t%s' % a)
            x = len(line.alleles)-1
            while x < numAltAlleles:
                outfile.write('\t')
                x += 1
        else:
            outfile.write('\t%s' % ','.join(line.alleles[1:]))
        
        line.extractQual()
        outfile.write("\t%f" % line.qual)
        
        line.extractFilters()
        if separateInfoFields:
            for f in line.filters:
                outfile.write('\t%s' % f)
            x = len(line.filters)
            while x < numFilters:
                outfile.write('\t')
                x += 1
        else:
            outfile.write('\t%s' % ','.join(line.filters))
        
        line.extractInfo()
        for i in infoOrder:
            if not line.info.has_key(i):
                if separateInfoFields:
                    for x in xrange(infoHeaders[i]):
                        outfile.write('\t')
                else:
                    outfile.write('\t')
            else:
                values = line.info[i]
                if not isinstance(values,list):
                    values = [values]
                for j,v in enumerate(values):
                    if v == None:
                        values[j] = i
                if separateInfoFields:
                    for v in values:
                        outfile.write('\t%s' % v)
                    x = len(values)
                    while x < infoHeaders[i]:
                        outfile.write('\t')
                        x += 1
                else:
                    outfile.write('\t%s' % ','.join(values))
        
        if includeGenotypes:
            line.extractFormat()
            line.extractGenotypes()
            for i,p in enumerate(peopleOrder):
                allele0,allele1,phased,attrs = line.genotypes[i]
                if allele0 == None:
                    allele0 = '.'
                elif not numberAlleles:
                    allele0 = line.alleles[allele0]
                if allele1 == None:
                    allele1 = '.'
                elif not numberAlleles:
                    allele1 = line.alleles[allele1]
                outfile.write('\t%s\t%s' % (allele0,allele1))
                if includeGenotypeAttributes:
                    outfile.write('\t%s' % ('Y' if phased else 'N'))
                    for j,f in enumerate(formatOrder):
                        if not f in line.format:
                            if separateInfoFields:
                                for x in xrange(formatHeaders[f]):
                                    outfile.write('\t')
                            else:
                                outfile.write('\t')
                        else:
                            attrIndex = line.format.index(f)-1
                            if attrIndex >= len(attrs):
                                if separateInfoFields:
                                    for x in xrange(formatHeaders[f]):
                                        outfile.write('\t')
                                else:
                                    outfile.write('\t')
                            else:
                                values = attrs[attrIndex].split(',')
                                if separateInfoFields:
                                    for v in values:
                                        outfile.write('\t%s' % v)
                                    x = len(values)
                                    while x < formatHeaders[f]:
                                        outfile.write('\t')
                                        x += 1
                                else:
                                    outfile.write('\t%s' % ','.join(values))
        outfile.write("\n")
    infile.close()
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a .csv file from the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO fields of a .vcf file. Also optionally includes columns for genotypes, or even their attributes.')
    parser.add_argument('--in', type=str, dest="infile", required=True,
                        help='Path to .vcf file')
    parser.add_argument('--out', type=str, dest="outfile", required=True,
                        help='Path to .csv file')
    parser.add_argument('--include_genotypes', type=str, dest="include_genotypes", nargs="?", const="True", default="False",
                        help='Include two columns for every individual, which will contain their called genotypes.')
    parser.add_argument('--include_genotype_attributes', type=str, dest="include_genotype_attributes", nargs="?", const="True", default="False",
                        help='Overrides --include_genotypes, and includes both genotypes and any values associated with them. This will probably be a very big file.')
    parser.add_argument('--separate_info_fields', type=str, dest="separate_info_fields", nargs="?", const="True", default="False",
                        help='When an INFO field has multiple comma-delimited values, split them into separate columns.')
    parser.add_argument('--numbered_alleles', type=str, dest="numbered_alleles", nargs="?", const="True", default="False",
                        help='Use the original numbered alleles .vcf-style instead of letters.')
    parser.add_argument('--ignore', type=str, dest="ignore_fields", nargs="+",
                        help='Explicitly remove specific columns from the output.')
    
    args = parser.parse_args()
    run(args)