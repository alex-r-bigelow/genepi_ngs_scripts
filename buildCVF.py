#!/usr/bin/env python
import argparse, os, datetime, math
from index_kgp import vcfLine, kgpInterface, countingDict

class infoDetails:
    def __init__(self, id, maxCategories, separateInfoFields):
        self.ranges = []    # nth column : (low,high) or None, indicating that the column exists, but there are no numerical values
        self.categories = []    # nth column : set(possible keys) or None, indicating that the column exists, but there are no strings
        
        self.globalRange = None
        self.globalCategories = None
        
        self.id = id
        self.maxCategories = maxCategories
        self.separateInfoFields = separateInfoFields
        
        self.maxedOut = False
    
    def addArbitraryValue(self, fields):
        if not isinstance(fields,list):
            fields = [fields]
        for i,f in enumerate(fields):
            try:
                f = float(f)
                if math.isnan(f):
                    self.addCategory('NaN',i)
                elif math.isinf(f):
                    self.addCategory("Inf", i)
                else:
                    self.addValue(f, i)
            except ValueError:
                self.addCategory(f, i)
    
    def addValue(self, v, i):
        while i >= len(self.ranges):
            self.ranges.append(None)
        if self.ranges[i] == None:
            self.ranges[i] = (v,v)
        else:
            self.ranges[i] = (min(v,self.ranges[i][0]),max(v,self.ranges[i][1]))
        if self.globalRange == None:
            self.globalRange = (v,v)
        else:
            self.globalRange = (min(v,self.globalRange[0]),max(v,self.globalRange[1]))
    
    def addCategory(self, v, i):
        while i >= len(self.categories):
            self.categories.append(None)
        if self.categories[i] == None:
            self.categories[i] = set()
        if self.separateInfoFields:
            if len(self.categories[i]) > self.maxCategories:
                self.maxedOut = True
                return
        else:
            if len(self.globalCategories) > self.maxCategories:
                self.maxedOut = True
                return
        self.categories[i].add(v)
        if self.globalCategories == None:
            self.globalCategories = set()
        self.globalCategories.add(v)
    
    def hasCategory(self, value, i=None):
        if i == None or not self.separateInfoFields:
            return value in self.globalCategories
        else:
            return value in self.categories[i]
    
    def getPragmas(self):
        results = []
        for i in xrange(max(len(self.ranges),len(self.categories))):
            pragmaString = "#\t%s %i\t" % (self.id, i+1)
            if self.maxedOut:
                pragmaString += "IGNORE"
            else:
                if self.separateInfoFields:
                    if self.ranges[i] == None:
                        if self.categories[i] == None:
                            continue
                        else:
                            pragmaString += "CATEGORICAL\t" + "\t".join(sorted(self.categories[i]))
                    else:
                        if self.categories[i] == None:
                            pragmaString += "NUMERIC\t%f\t%f" % self.ranges[i]
                        else:
                            pragmaString += "MIXED\t%f\t%f\t" % self.ranges[i] + "\t".join(sorted(self.categories[i]))
                else:
                    if self.globalRange == None:
                        if self.globalCategories == None:
                            continue
                        else:
                            pragmaString += "CATEGORICAL\t" + "\t".join(sorted(self.globalCategories))
                    else:
                        if self.globalCategories == None:
                            pragmaString += "NUMERIC\t%f\t%f" % self.globalRange
                        else:
                            pragmaString += "MIXED\t%f\t%f\t" + self.ranges[i] + "\t".join(sorted(self.categories[i]))
            results.append(pragmaString)
        return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a .cvf file from the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO fields of a .vcf file')
    parser.add_argument('--in', type=str, dest="infile",
                        help='Path to .vcf file')
    parser.add_argument('--out', type=str, dest="outfiles", nargs="+",
                        help='Path to .cvf file')
    parser.add_argument('--max_strings', type=int, dest="max_strings", nargs="?", const=50, default=50,
                        help='Maximum number of strings a categorical INFO field can have before it\'s automatically marked as IGNORE')
    parser.add_argument('--separate_info_fields', type=str, dest="separate_info_fields", nargs="?", const="False", default="False",
                        help='When an INFO field has multiple comma-delimited values, count possible ranges and categories separately')
    
    args = parser.parse_args()
    
    separateInfoFields = args.separate_info_fields.strip().lower() == "true"
    
    maxAlleles = 0
    maxFilters = 0
    pragmaFilterTags = set()
    filterColumn = infoDetails("FILTER", args.max_strings, separateInfoFields)
    maxQual = 0.0
    infoFields = {}
    # TODO: get the numeric ranges, all valid categorical values
    
    print "Counting columns, ranges, values...",
    infile = open(args.infile,'r')
    for line in infile:
        line = line.strip()
        if len(columns) <= 1:
            continue
        elif line.startswith("#"):
            if line.startswith("##INFO"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                if infoFields.has_key(newTag):
                    raise Exception("Duplicate INFO ID or use of reserved ID:\t%s" % newTag)
                infoFields[newTag] = infoDetails(newTag, args.max_strings,separateInfoFields)
            elif line.startswith("##FILTER"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                if newTag in pragmaFilterTags:
                    raise Exception("Duplicate FILTER ID: %s" % newTag)
                pragmaFilterTags.add(newTag)
                if not separateInfoFields:
                    filterColumn.addCategory(newTag, 0)
            continue
        else:
            line = vcfLine(line.split('\t'))
            line.extractAlleles()
            line.extractQual()
            line.extractFilters()
            line.extractInfo()
            maxAlleles = max(maxAlleles,len(line.alleles))
            maxQual = max(maxQual,line.qual)
            maxFilters = max(maxFilters,len(line.filters))
            for i,f in enumerate(line.filters):
                if not f in pragmaFilterTags:
                    raise Exception("Missing ##FILTER pragma for: %s" % f)
                filterColumn.addCategory(f,i)
            
            for k,v in line.info.iteritems():
                if k == 'ALT' or k == 'FILTER':
                    continue
                elif not infoFields.has_key(k):
                    raise Exception("Missing ##INFO pragma for: %s" % k)
                else:
                    infoFields[k].addArbitraryValue(v)
    infile.close()
    
    print "Creating file..."
    outfile = open(args.outfile, 'w')
    
    outfile.write("##\t%s created from %s on %s\n" % (args.outfile,args.str(infile,datetime.datetime.now())))
    outfile.write("#\tChromosome\tCHR\n")
    outfile.write("#\tPosition\tPOS\n")
    outfile.write("#\tID\tID\n")
    for i in xrange(len(numColumns['ALT'])):
        outfile.write("#\tAllele %i\tIGNORE\n" % (i+1))
    outfile.write("#\tQUAL\tNUMERIC\t0.0\t%f\n" % maxQual)
    for k,v in sorted(numColumns.iteritems()):
        # TODO: continue here!
    
    outfile.write('Chromosome\tPosition')
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i'%i)
    for i in xrange(numAlleles):
        outfile.write('\t# Sharing Allele %i'%i)
    outfile.write('\t# Samples w/data')
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i Frequency'%i)
    for i in xrange(numAlleles):
        outfile.write('\tAllele %i KGP Frequency'%i)
    
    total = 0
    noData = 0
    passedKGP = 0
    sharingCounts = countingDict()
    sharingFractions = {}
        
    for line,kgpLine in kgp.iterateVcf(args.infile,tickFunction=tick,numTicks=1000):
        outfile.write('\n')
        
        # write chromosome, position
        outfile.write('%s\t%i'%(line.chromosome,line.position))
        
        line.extractAlleles()
        
        for a in line.alleles:
            outfile.write('\t%s'% a)
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t-')
        
        # write sharing per allele
        sharing = line.getSharing()
        for a in line.alleles:
            outfile.write('\t%i'%sharing[a])
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t-')
        
        # write sharing denominator
        numberWithData = line.numberOfSamplesWithData()
        outfile.write('\t%i'%numberWithData)
        
        # write allele frequencies in file
        frequencies = line.getAlleleFrequencies()
        for a in line.alleles:
            outfile.write('\t%f'%frequencies[a])
        for x in xrange(numAlleles-len(line.alleles)):
            outfile.write('\t-')
        
        # write kgp allele frequencies
        kgpFrequencies = None
        if kgpLine == None:
            for a in line.alleles:
                outfile.write('\tND')
            for x in xrange(numAlleles-len(line.alleles)):
                outfile.write('\t-')
        else:
            kgpLine.extractGenotypes()
            kgpFrequencies = kgpLine.getAlleleFrequencies()
            for a in line.alleles:
                if kgpFrequencies.has_key(a):
                    outfile.write('\t%f'%kgpFrequencies[a])
                else:
                    outfile.write('\tND')
            for x in xrange(numAlleles-len(line.alleles)):
                outfile.write('\t-')
        
        # count stuff for final output
        total += 1
        if kgpLine == None or numberWithData == 0:
            noData += 1
        else:
            lowestFreq = 1.0
            minorAllele = None
            for a in line.alleles:
                if kgpFrequencies.has_key(a):
                    if kgpFrequencies[a] < lowestFreq:
                        lowestFreq = kgpFrequencies[a]
                        minorAllele = a
            if lowestFreq <= 0.01:
                passedKGP += 1
                
                proportion = sharing[minorAllele]/float(numberWithData)
                sharingCounts[proportion] += 1
                sharingFractions[proportion] = (sharing[minorAllele],numberWithData)
    
    infile.close()
    outfile.close()
    print ""
    print "Done"
    print ""
    print "Total number of variants: %i" % total
    print "# of variants with at least some data from both sources: %i" % (total-noData)
    print "# left if variants with no allele <= 0.01 are removed: %i" % passedKGP
    temptotal = 0
    for p,c in sorted(sharingCounts.iteritems()):
        print "# left if variants sharing < %i/%i of the minor KGP allele are removed: %i" % (sharingFractions[p][0],sharingFractions[p][1],passedKGP-temptotal)
        temptotal += c