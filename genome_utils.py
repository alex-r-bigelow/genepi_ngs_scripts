#!/usr/bin/env python
import os, gzip, math

MAX_INFO_STRINGS=40

chromosomeOrder = ['chr1',
                   'chr2',
                   'chr3',
                   'chr4',
                   'chr5',
                   'chr6',
                   'chr7',
                   'chr8',
                   'chr9',
                   'chr10',
                   'chr11',
                   'chr12',
                   'chr13',
                   'chr14',
                   'chr15',
                   'chr16',
                   'chr17',
                   'chr18',
                   'chr19',
                   'chr20',
                   'chr21',
                   'chr22',
                   'chrX',
                   'chrY',
                   'chrM']
chromosomesToKeep = set(chromosomeOrder)
chromosomeRank = {} # for comparing variants, this is a little faster
for i,c in enumerate(chromosomeOrder):
    chromosomeRank[c] = i

class genomeException(Exception):
    pass

def parsePopulations(path):
    with open(path,'rb') as infile:
        populations = {}
        header = None
        for line in infile:
            if line.startswith('#\t'):
                if header == None:
                    header = line.strip().split('\t')[1:]
                else:
                    columns = line.strip().split('\t')[1:]
                    for i,c in enumerate(columns):
                        if c != "":
                            if not populations.has_key(header[i]):
                                populations[header[i]] = []
                            populations[header[i]].append(c)
        return (populations,header)

def standardizeChromosome(chrom):
    if not chrom.lower().startswith("chr"):
        return "chr" + chrom
    else:
        return "chr" + chrom[3:]

class countingDict(dict):
    def __missing__(self, key):
        returnValue = 0
        self[key] = returnValue
        return returnValue

class vcfLine:
    @staticmethod
    def constructLine(chromosome,position,name=".",alleles=["N","N"],info={},qual=0.0,filters=["."],format=["GT"],number_of_genotypes=0):
        temp = vcfLine([])
        
        chromosome = standardizeChromosome(chromosome)
        temp.chromosome = chr
        temp.columns.append(chr)
        
        temp.position = position
        temp.columns.append(str(position))
        
        temp.name = name
        temp.columns.append(name)
        
        temp.alleles = alleles
        temp.columns.append(alleles[0])
        temp.columns.append(",".join(alleles[1:]))
        
        temp.qual = qual
        temp.columns.append(str(qual))
        
        temp.filters = filters
        temp.columns.append(";".join(sorted(filters)))
        
        temp.info = info
        infostrs = []
        for k,v in info.iteritems():
            if v == None:
                infostrs.append(k)
            elif isinstance(v,list):
                infostrs.append(k + "=" + ','.join(v))
            else:
                infostrs.append(k + "=" + v)
        temp.columns.append(";".join(sorted(infostrs)))
        
        temp.format = format
        temp.columns.append(":".join(format))
        
        temp.genotypes = {}
        for i in xrange(number_of_genotypes):
            temp.genotypes[i] = (None,None)
            temp.columns.append("./.")
        
        return temp
    
    def __init__(self, columns):
        self.columns = columns
    
    def extractChrAndPos(self):
        if not hasattr(self,'chromosome'):
            self.chromosome = standardizeChromosome(self.columns[0])
            self.position = int(self.columns[1])
            self.name = self.columns[2]
    
    def extractAlleles(self):
        if not hasattr(self, 'alleles'):
            self.alleles = self.columns[4].split(',')
            self.alleles.insert(0,self.columns[3])
    
    def extractInfo(self):
        if not hasattr(self, 'info'):
            self.info = {}
            fields = self.columns[7].split(';')
            for i in fields:
                if '=' in i:
                    values = i.split('=')
                    key = values[0]
                    value = values[1]
                    if ',' in value:
                        value = value.split(',')
                elif i == "":
                    continue
                else:
                    key = i
                    value = None
                
                self.info[key] = value
    
    def extractQual(self):
        if not hasattr(self, 'qual'):
            self.qual = float(self.columns[5])
    
    def extractFilters(self):
        if not hasattr(self, 'filters'):
            self.filters = self.columns[6].split(';')
    
    def extractFormat(self):
        if not hasattr(self, 'format'):
            self.format = self.columns[8].split(':')
    
    def extractGenotypes(self, indices=None):
        if indices == None:
            indices = range(len(self.columns[9:]))
        if not hasattr(self, 'genotypes'):
            self.genotypes = {}
        for i in indices:
            if not self.genotypes.has_key(i):
                genotype = self.columns[i+9].split(':')
                if "|" in genotype[0]:
                    phased = True
                    alleles = genotype[0].split("|")
                else:
                    phased = False
                    alleles = genotype[0].split("/")
                if alleles[0] == ".":
                    allele0 = None
                else:
                    allele0 = int(alleles[0])
                if len(alleles) == 1:
                    allele1 = -1
                elif alleles[1] == ".":
                    allele1 = None
                else:
                    allele1 = int(alleles[1])
                
                if allele0 == None or allele1 == None:
                    assert allele0 == None and allele1 == None
                
                self.genotypes[i] = (allele0,allele1,phased,genotype[1:])
    
    def reorderAlleles(self, alleleScores, highToLow=True):
        ''' Reorder the alleles (and genotype numbers) based on some positive score for each allele (usually a background allele frequency).
        At the moment I'm too lazy to implement reordering of genotype attributes as well (it gets a little hairy with specific
        fields), so if you use this method, all genotype attributes will be removed from the output .vcf file.
        Note as well that "REF" and "ALT" may be a little misleading as the "REF" allele will no longer be the reference
        allele (but ALT would become the true minor allele in the case of allele frequencies and highToLow sorting). Alleles
        not in alleleScores are assumed to have a score of 0 '''
        self.extractAlleles()
        
        scoremap = sorted(alleleScores.iteritems(),key=lambda x:x[1],reverse=highToLow)
        alleleMap = {}
        newAlleles = []
        
        for i,(a,s) in enumerate(scoremap):
            if a in self.alleles:
                alleleMap[self.alleles.index(a)] = i
                newAlleles.append(a)
        
        for i,a in enumerate(self.alleles):
            if a not in newAlleles:
                alleleMap[i] = len(newAlleles)
                newAlleles.append(a)
        
        self.alleles = newAlleles
        
        self.extractGenotypes()
        for i,g in self.genotypes.iteritems():
            allele0,allele1,phased,genAttrs = g
            if allele0 != None:
                allele0 = alleleMap[allele0]
            if allele1 != None:
                allele1 = alleleMap[allele1]
            
            self.genotypes[i] = (allele0,allele1,phased,[])
    
    def __repr__(self):
        outline = ""
        
        if hasattr(self,'chromosome'):
            outline += self.chromosome + '\t'
        else:
            outline += standardizeChromosome(self.columns[0]) + '\t'
        
        if hasattr(self,'position'):
            outline += '%i\t' % self.position
        else:
            outline += self.columns[1] + '\t'
        
        if hasattr(self,'name'):
            outline += self.name + '\t'
        else:
            outline += self.columns[2] + '\t'
        
        if hasattr(self,'alleles'):
            outline += self.alleles[0] + '\t' + ','.join(self.alleles[1:]) + '\t'
        else:
            outline += self.columns[3] + '\t' + self.columns[4] + '\t'
        
        if hasattr(self,'qual'):
            outline += '%f\t' % self.qual
        else:
            outline += self.columns[5] + '\t'
        
        if hasattr(self,'filters'):
            outline += ";".join(self.filters) + '\t'
        else:
            outline += self.columns[6] + '\t'
        
        if hasattr(self,'info'):
            infostrs = []
            for k,v in self.info.iteritems():
                if v == None:
                    infostrs.append(k)
                elif isinstance(v,list):
                    infostrs.append(k + "=" + ','.join(v))
                else:
                    infostrs.append(k + "=" + v)
            outline += ";".join(sorted(infostrs)) + '\t'
        else:
            outline += self.columns[7] + '\t'
        
        if hasattr(self,'format'):
            outline += ":".join(self.format) + '\t'
        else:
            outline += self.columns[8] + '\t'
        
        for i in xrange(len(self.columns[9:])):
            if hasattr(self,'genotypes') and self.genotypes.has_key(i):
                if self.genotypes[i][0] == None:
                    outline += "."
                else:
                    outline += str(self.genotypes[i][0])
                if self.genotypes[i][2]:
                    outline += "|"
                else:
                    outline += "/"
                if self.genotypes[i][1] == None:
                    outline += "."
                elif not self.genotypes[i][1] == -1:
                    outline += str(self.genotypes[i][1])
                if len(self.genotypes[i][3]) > 0:
                    outline += ":" + ":".join(self.genotypes[i][3])
                outline += "\t"
            else:
                outline += self.columns[i+9] + "\t"
        
        return outline + '\n'

class infoDetails:
    def __init__(self, id, maxCategories, countSeparate):
        self.ranges = []    # nth column : (low,high) or None, indicating that the column exists, but there are no numerical values
        self.categories = []    # nth column : set(possible keys) or None, indicating that the column exists, but there are no strings
                
        self.name = id
        self.maxCategories = maxCategories
        self.countSeparate = countSeparate
        
        self.numColumns = 1
        self.maxedOut = False
    
    def addArbitraryValue(self, fields):
        if not isinstance(fields,list):
            fields = [fields]
        self.numColumns = max(len(fields),self.numColumns)
        for i,f in enumerate(fields):
            if f == None:
                continue
            try:
                f = float(f)
                if math.isnan(f):
                    self.addCategory('NaN',i)
                elif math.isinf(f):
                    self.addCategory("Inf",i)
                else:
                    self.addValue(f,i)
            except ValueError:
                self.addCategory(f,i)
    
    def addValue(self, v, i):
        if not self.countSeparate:
            i = 0
        while i >= len(self.ranges):
            self.ranges.append(None)
        
        if self.ranges[i] == None:
            self.ranges[i] = (v,v)
        else:
            self.ranges[i] = (min(v,self.ranges[i][0]),max(v,self.ranges[i][1]))
    
    def addCategory(self, v, i):
        if self.maxedOut:
            return
        
        if not self.countSeparate:
            i = 0
        while i >= len(self.categories):
            self.categories.append(None)
        
        if self.categories[i] == None:
            self.categories[i] = set()
        self.categories[i].add(v)
        
        if self.maxCategories > 0 and len(self.categories[i]) > self.maxCategories:
            self.maxedOut = True
    
    def hasCategory(self, value, i=0):
        return value in self.categories[i]
    
    def getPragmas(self, override=None):
        results = []
        for i in xrange(self.numColumns):
            if i == 0:
                pragmaString = "#\t%s\t" % self.name
            else:
                pragmaString = "#\t%s %i\t" % (self.name, i+1)
            
            if override != None:
                pragmaString += override
            elif self.maxedOut:
                pragmaString += "IGNORE"
                results.append(pragmaString)
                continue
            else:
                if self.ranges[i] == None:
                    if self.categories[i] == None:
                        pragmaString += "IGNORE"
                    else:
                        pragmaString += "CATEGORICAL"
                else:
                    if self.categories[i] == None:
                        pragmaString += "NUMERIC"
                    else:
                        pragmaString += "MIXED"
            
            if self.ranges[i] != None:
                # TODO: allow rounding to nearest tens...
                pragmaString += "\t(%f,%f)" % self.ranges[i]
            
            if self.categories[i] != None:
                # TODO: sort categories by their frequency?
                pragmaString += "\t" + "\t".join(sorted(self.categories[i]))
            
            results.append(pragmaString)
        return results

class kgpInterface:
    BYTES_TO_ITERATE=8*4096
    def __init__(self, dataPath, popPath):
        '''
        Creates an interface to 1000 genomes .vcf.gz files; these should be downloaded
        in a single directory (dirpath should be the path to that directory). This interface utilizes
        the default file names and sort ordering; don't modify what you download!
        
        You should also have a KGP_populations.txt file around somewhere describing the
        population structure in the 1000 genomes (popPath) - see calcStats.py --help for more
        details about its format
        '''
        self.populations = parsePopulations(popPath)[0]
        self.populationIndices = {}
        self.files = {}
        if dataPath != None:
            self.valid = True
            for dirname, dirnames, filenames in os.walk(dataPath):
                for filename in filenames:
                    if filename.endswith('.vcf.gz') and 'chr' in filename:
                        chrname = filename[filename.find("chr"):]
                        chrname = chrname[:chrname.find(".")]
                        fullpath = os.path.join(dirname,filename)
                        assert chrname in chromosomeOrder
                        self.files[chrname]=gzip.open(fullpath,'rb')
            for line in self.files.itervalues().next(): # just grab one of the files
                if line.startswith("#CHROM"):
                    self.header = line.strip().split('\t')
                    for p,individuals in self.populations.iteritems():
                        self.populationIndices[p] = []
                        for i in individuals:
                            self.populationIndices[p].append(self.header.index(i)-9)
                    break
            self.startAtZero()
        else:
            self.valid = False
    
    def startAtZero(self):
        for f in self.files.itervalues():
            f.seek(0)
    
    def iterate(self):
        self.startAtZero()
        return self._iterate()
    
    def _iterate(self):
        for c in chromosomeOrder:
            if not self.files.has_key(c):
                continue
            people = None
            for line in self.files[c]:
                if line.startswith('#CHROM'):
                    people = line.strip().split('\t')[9:]
                    continue
                elif people == None:
                    continue
                else:
                    yield(line,people)
    
    def iterateVcf(self, vcfPath, tickFunction=None, numTicks=100):
        ''' Useful for iterating through a sorted .vcf file and finding matches in KGP; the vcf file should be
        base pair position-ordered (the chromosome order is irrelevant) '''
        infile = open(vcfPath,'rb')
        tickInterval = os.path.getsize(vcfPath)/numTicks
        # We take advantage of the fact that the KGP .vcf files are bp-ordered
        self.startAtZero()
        return self._iterateVcf(infile, tickFunction, tickInterval)            
    
    def _iterateVcf(self, infile, tickFunction, tickInterval):
        nextTick = 0
        kgpLines = {}
        for f in self.files.iterkeys():
            kgpLines[f] = None
                
        for vline in infile:
            vline = vline.strip()
            if len(vline) <= 1 or vline.startswith("#"):
                continue
            if tickFunction != None and infile.tell() >= nextTick:
                nextTick += tickInterval
                tickFunction()
            vline = vcfLine(vline.split('\t'))
            vline.extractChrAndPos()
            
            # If we're missing data for a particular chromosome (e.g. chrMT, etc), just harmlessly return that that line is missing
            if not self.valid or not kgpLines.has_key(vline.chromosome):
                yield (vline,None)
                continue
            
            # continue through the KGP file (we assume the .vcf file is sorted by position, chromosome doesn't matter)
            # until we match or pass the .vcf line
            eof = False
            while kgpLines[vline.chromosome] == None or kgpLines[vline.chromosome].position < vline.position:
                text = self.files[vline.chromosome].readline()
                if text == '':
                    eof = True
                    break
                text = text.strip()
                if len(text) <= 1 or text.startswith('#'):
                    continue
                else:
                    kgpLines[vline.chromosome] = vcfLine(text.split('\t'))
                    kgpLines[vline.chromosome].extractChrAndPos()
                    assert kgpLines[vline.chromosome].chromosome == vline.chromosome
            if eof:
                yield (vline,None)
                continue
            elif kgpLines[vline.chromosome].position == vline.position:
                yield (vline,kgpLines[vline.chromosome])
                continue
            else:
                yield (vline,None)
                continue
        # way out here, the .vcf file is depleted; close the file and we're done
        infile.close()
        raise StopIteration

class bedLine:
    def __init__(self, columns):
        self.columns = columns
        self.chromosome = standardizeChromosome(self.columns[0])
        self.start = int(self.columns[1])+1 # Internally, I use 1-based coordinates like .vcf files do
        self.stop = int(self.columns[2])+1
        if len(self.columns) > 3:
            self.name = self.columns[3]
        if len(self.columns) > 4:
            self.score = float(self.columns[4])
    
    def contains(self, position):
        return position >= self.start and position < self.stop
    
    def __repr__(self):
        outline = "%s\t%i\t%i" % (self.chromosome,self.start-1,self.stop-1)
        
        if hasattr(self.name):
            outline += "\t" + self.name
        else:
            outline += "\t" + self.columns[3]
        
        if hasattr(self.name):
            outline += "\t" + str(self.score)
        else:
            outline += "\t" + self.columns[4]
        
        if len(self.columns) > 5:
            outline += "\t" + "\t".join(self.columns[5:])
        
        return outline + "\n"