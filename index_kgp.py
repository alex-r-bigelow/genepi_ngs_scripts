#!/usr/bin/env python
import pickle, argparse, sys, os

class countingDict(dict):
    def __missing__(self, key):
        returnValue = 0
        self[key] = returnValue
        return returnValue

class vcfLine:
    def __init__(self, columns):
        self.columns = columns
    
    def extractChrAndPos(self):
        if not hasattr(self,'chromosome'):
            self.chromosome = self.columns[0].lower()
            if not self.chromosome.startswith('chr'):
                self.chromosome = 'chr' + self.chromosome
            self.position = int(self.columns[1])
            self.id = self.columns[2]
    
    def extractAlleles(self):
        if not hasattr(self, 'alleles'):
            self.alleles = self.columns[4].split(',')
            self.alleles.insert(0,self.columns[3])
    
    def extractInfo(self):
        if not hasattr(self, 'info'):
            self.info = {}
            for i in self.columns[7].split(';'):
                if '=' in i:
                    values = i.split('=')
                    key = i[0]
                    value = i[1]
                    if ',' in value:
                        value = value.split(',')
                else:
                    key = i
                    value = None
                
                self.info[key] = value
    
    def extractQual(self):
        if not hasattr(self, 'qual'):
            self.qual = float(columns[5])
    
    def extractFilters(self):
        if not hasattr(self, 'filters'):
            self.filters = columns[6].split(';')
    
    def extractFormat(self):
        if not hasattr(self, 'format'):
            self.format = self.columns[8].split(':')
    
    def extractGenotypes(self, indices=None):
        if indices == None:
            indices = range(len(self.columns))[9:]
        if not hasattr(self, 'genotypes'):
            self.genotypes = {}
        for i in indices:
            if not self.genotypes.has_key(i):
                genotype = self.columns[i].split(':')
                if "|" in genotype[0]:
                    alleles = genotype[0].split("|")
                else:
                    alleles = genotype[0].split("/")
                if alleles[0] == ".":
                    allele0 = None
                else:
                    allele0 = int(alleles[0])
                if alleles[1] == ".":
                    allele1 = None
                else:
                    allele1 = int(alleles[1])
                
                if allele0 == None or allele1 == None:
                    assert allele0 == None and allele1 == None
                
                self.genotypes[i] = (allele0,allele1)
    
    def numberOfSamplesWithData(self, indices=None):
        self.extractGenotypes(indices)
        if indices == None:
            indices = self.genotypes.iterkeys()
        count = 0
        for i in indices:
            if self.genotypes[i][0] != None:
                count += 1
        return count
    
    def getAlleleFrequencies(self, indices=None):
        '''
        population should be either a string that is a key in self.populations or a set or list of integer column numbers with which to query line
        '''
        self.extractAlleles()
        self.extractGenotypes(indices)
        if indices == None:
            indices = self.genotypes.iterkeys()
        
        count = 0.0
        matches = countingDict()
        for i in indices:
            allele0, allele1 = self.genotypes[i]
            if allele0 != None:
                count += 2.0
                matches[allele0] += 1
                matches[allele1] += 1
        
        if count == 0.0:
            return dict([(a,float('Inf')) for a in self.alleles])
        else:
            return dict([(a,matches.get(i,0)/count) for i,a in enumerate(self.alleles)])
    
    def getSharing(self, indices=None):
        '''
        population should be either a string that is a key in self.populations or a set or list of integer column numbers with which to query vcfLine
        '''
        
        self.extractAlleles()
        self.extractGenotypes(indices)
        if indices == None:
            indices = self.genotypes.iterkeys()
        
        counts = countingDict()
        for i in indices:
            allele0,allele1 = self.genotypes[i]
            
            if allele0 != None:
                if allele0 == allele1:
                    counts[allele0] += 1
                else:
                    counts[allele0] += 1
                    counts[allele1] += 1
        return dict([(a,counts.get(i,0)) for i,a in enumerate(self.alleles)])
    
    def __str__(self):
        outline = ""
        
        if hasattr(self,'chromosome'):
            outline += self.chromosome + '\t'
        else:
            outline += self.columns[0] + '\t'
        
        if hasattr(self,'position'):
            outline += '%i\t' % self.position
        else:
            outline += self.columns[1] + '\t'
        
        if hasattr(self,'id'):
            outline += self.id + '\t'
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
        
        if len(columns) >= 9:
            # TODO: be fancier with the genotypes...?
            outline += '\t'.join(self.columns[9:])
        
        return outline + '\n'

class kgpInterface:
    BYTES_TO_ITERATE=8*4096
    def __init__(self, dirpath):
        '''
        Creates an interface to 1000 genomes .vcf files; these should be downloaded and decompressed
        in a single directory (dirpath should be the path to that directory). This interface utilizes
        the default file names and sort ordering; don't modify what you download!
        
        You should also have run this app at least once to build KGP.index in that directory before
        using this class.
        '''
        infile = open(os.path.join(dirpath,"KGP.index"),'rb')
        self.populations = pickle.load(infile)
        temp = pickle.load(infile)
        self.fileSizes = pickle.load(infile)
        self.files = {}
        for k,v in temp.iteritems():
            self.files[k]=open(v,'r')
        infile.close()
    
    def iterateVcf(self, vcfPath, tickFunction=None, numTicks=1000):
        ''' Useful for iterating through a sorted .vcf file and finding matches in KGP; the vcf file should be
        base pair position-ordered (the chromosome order is irrelevant) '''
                
        # We take advantage of the fact that the KGP .vcf files are bp-ordered
        for f in self.files.itervalues():
            f.seek(0)
        
        infile = open(vcfPath,'r')
        
        tickInterval = os.path.getsize(vcfPath)/numTicks
        
        return self._iterateVcf(infile, tickFunction, tickInterval)
    
    def _iterateVcf(self, infile, tickFunction, tickInterval):
        nextTick = 0
        lastKGPline = None
                
        for vline in infile:
            vline = vline.strip()
            if len(vline) <= 1 or vline.startswith("#"):
                continue
            if tickFunction != None and infile.tell() >= nextTick:
                nextTick += tickInterval
                tickFunction()
            vline = vcfLine(vline.split('\t'))
            vline.extractChrAndPos()
            
            if not self.files.has_key(vline.chromosome):
                yield (vline,None)
                continue
            
            # special case: if the last legit comparison was no match, we will have a line left over from KGP that could still
            # potentially match
            if lastKGPline != None:
                if vline.chromosome != lastKGPline.chromosome:
                    # wierd little corner case at the end of a chromosome; the last pass could be irrelevant
                    lastKGPline = None
                elif vline.position < lastKGPline.position:
                    # Another no match; leave lastKGPline alone
                    yield (vline,None)
                    continue
                elif vline.position == lastKGPline.position:
                    # A match! Make sure to use up lastKGPline
                    yield (vline,lastKGPline)
                    lastKGPline = None
                    continue
                else:
                    # okay the kgp line won't ever be matched, so we can continue with iterating
                    lastKGPline = None
            
            # keep reading lines from KGP from where we left off last until we find/pass a match
            yieldedSomething = False
            for line in self.files[vline.chromosome]:
                line = line.strip()
                if len(line) <= 1 or line.startswith("#"):
                    continue
                line = vcfLine(line.split('\t'))
                line.extractChrAndPos()
                if vline.position < line.position:
                    # this line is too early; there's no match, but we just read a line from KGP so we need to store the special case
                    lastKGPline = line
                    yieldedSomething = True
                    yield (vline,None)
                elif vline.position == line.position:
                    # a match!
                    yieldedSomething = True
                    yield (vline,line)
                else:
                    # this is a little funky; we want to keep going if we haven't found/passed a match in KGP
                    lastByte = self.files[vline.chromosome].tell()
                    continue
                # ... but if we actually yielded something, we want to move on to the next .vcf line
                break
            # out here, we've reached the end of KGP; yield no match until the .vcf file is depleted
            if not yieldedSomething:
                yield (vline,None)
        # way out here, the .vcf file is depleted; close the file and we're done
        infile.close()
        raise StopIteration
    
    def getVcfLines(self, chromosome, start, stop):
        # Find a reasonable range of the file to iterate through by first
        # using our hash table
        lowByte = self.fileSizes[chromosome][0]
        highByte = self.fileSizes[chromosome][1]
        # now binary search until the space is iteratable
        while highByte-lowByte > kgpInterface.BYTES_TO_ITERATE:
            self.files[chromosome].seek((lowByte+highByte)/2)
            dummy = self.files[chromosome].readline()   # throw away what's left of this line
            nextByte = self.files[chromosome].tell()
            dummy = self.files[chromosome].readline()
            columns = dummy.strip().split('\t')
            pos = int(columns[1])
            if pos >= start:
                highByte = nextByte
            else:
                lowByte = nextByte
        # now iterate until we hit stop or EOF, collecting lines once we pass start
        self.files[chromosome].seek(lowByte)
        linesToReturn = []
        for line in self.files[chromosome]:
            columns = line.strip().split('\t')
            pos = int(columns[1])
            if pos >= start:
                if pos >= stop:
                    return linesToReturn
                else:
                    linesToReturn.append(vcfLine(columns))
        return linesToReturn

if __name__ == "__main__":
    # If we're actually running this as a program, we want to build the index
    parser = argparse.ArgumentParser(description='Builds an index into the thousand genomes data.')
    parser.add_argument('--data', type=str, dest="data",
                        help='Directory housing decompressed .vcf files 1-22,X')
    parser.add_argument('--populations', type=str, dest="pop",
                        help="Directory containing .txt lists of KGP sample IDs; the file name should correspond to the population.")
    
    args = parser.parse_args()
    
    # {population name:set(sample ID)}
    populations = {}
    reversePopulations = {}
    
    for dirname, dirnames, filenames in os.walk(args.pop):
        for filename in filenames:
            if filename.endswith('.txt'):
                popName = filename.strip()[:-4]
                populations[popName] = {}
                infile = open(os.path.join(dirname,filename),'r')
                for line in infile:
                    populations[popName][line.strip()] = None
                    reversePopulations[line.strip()]=popName
                infile.close()
    
    # {chrN:file name}
    files = {}
    
    # {chrN:(first byte of real data,last byte in file)}
    fileSizes = {}
    
    for dirname, dirnames, filenames in os.walk(args.data):
        for filename in filenames:
            if filename.endswith('.vcf') and 'chr' in filename:
                chrname = filename[filename.find("chr"):]
                chrname = chrname[:chrname.find(".")]
                fullpath = os.path.join(dirname,filename)
                files[chrname]=fullpath
                
                infile = open(fullpath,'r')
                for line in infile:
                    if line.startswith("##") or len(line.strip()) <= 1:
                        continue
                    elif line.startswith("#"):
                        columns = line.strip().split('\t')
                        for i,c in enumerate(columns[9:]):
                            if reversePopulations.has_key(c):
                                if populations[reversePopulations[c]][c] == None:
                                    populations[reversePopulations[c]][c] = i+9
                                elif populations[reversePopulations[c]][c] != i+9:
                                    raise Exception("KGP files have %s in different columns: %i, %i" % (c,i+9,populations[reversePopulations[c]][c]))
                        fileSizes[chrname] = (infile.tell(),os.path.getsize(fullpath))
                        break
                    else:
                        raise Exception("KGP file %s has no header line"%filename)
                infile.close()
    
    filestr = os.path.join(args.data,"KGP.index")
    outfile = open(filestr,'wb')
    pickle.dump(populations,outfile)
    pickle.dump(files,outfile)
    pickle.dump(fileSizes,outfile)
    outfile.close()