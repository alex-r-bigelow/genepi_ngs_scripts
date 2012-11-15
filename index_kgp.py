#!/usr/bin/env python
import pickle, argparse, sys, os

class countingDict(dict):
    def __missing__(self, key):
        returnValue = 0
        self[key] = returnValue
        return returnValue

class vcfLine:
    def __init__(self, line):
        self.columns = line.strip().split('\t')
    
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
    
    def extractOther(self):
        if not hasattr(self, 'qual'):
            self.qual = float(columns[5])
            self.filters = columns[6].split(';')
            self.info = self.columns[7].split(';')
            self.format = self.columns[8].split(':')
    
    def extractGenotypes(self, indices=None):
        if indices == None:
            indices = range(len(self.columns))[9:]
        if not hasattr(self, 'genotypes'):
            self.genotypes = {}
        for i in self.indices:
            if not self.genotypes.has_key(i):
                genotype = self.columns[i].split(':')
                if "|" in genotype:
                    alleles = genotype.split("|")
                else:
                    alleles = genotype.split("/")
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
    
    def numberOfSamplesWithData(self, population):
        self.extractGenotypes(population)
        count = 0
        for i in population:
            if self.genotypes[i] != None:
                count += 1
        return count
    
    def getAlleleFrequencies(self, population):
        '''
        population should be either a string that is a key in self.populations or a set or list of integer column numbers with which to query line
        '''
        self.extractAlleles()
        self.extractGenotypes(population)
        
        count = 0.0
        matches = countingDict()
        for i in population:
            allele0, allele1 = self.genotypes[i]
            if allele0 != None:
                count += 2.0
                matches[allele0] += 1
                matches[allele1] += 1
        
        if count == 0.0:
            return dict([(a,float('Inf') for a in self.alleles)])
        else:
            return dict([(a,matches.get(i,0)/count) for i,a in enumerate(self.alleles)])
    
    def getSharing(self, population):
        '''
        population should be either a string that is a key in self.populations or a set or list of integer column numbers with which to query vcfLine
        '''
        
        self.extractAlleles()
        self.extractGenotypes(population)
        
        counts = countingDict()
        for i in population:
            allele0,allele1 = line.genotypes[i]
            
            if allele0 != None:
                if allele0 == allele1:
                    counts[allele0] += 1
                else:
                    counts[allele0] += 1
                    counts[allele1] += 1
        return dict([(a,counts.get(i,0)) for i,a in enumerate(self.alleles)])

class kgpInterface:
    def __init__(self, dirpath):
        infile = open(os.path.join(dirpath,"KGP.index"),'rb')
        self.level = pickle.load(infile)
        self.populations = pickle.load(infile)
        self.positions = pickle.load(infile)
        self.files = {}
        for k,v in pickle.load(infile):
            self.files[k]=open(v,'r')
        infile.close()
    
    def getVcfLines(self, chromosome, start, stop):
        self.files[chromosome].seek(self.positions[chromosome][int(start/self.level)])
        linesToReturn = []
        passedStart = False
        for line in self.files[chromosome]:
            columns = line.split('\t')
            pos = int(columns[1])
            if pos >= start:
                if pos >= stop:
                    return linesToReturn
                else:
                    linesToReturn.append(vcfLine(line))
        return linesToReturn

if __name__ == "__main__":
    # If we're actually running this as a program, we want to build the index
    parser = argparse.ArgumentParser(description='Builds an index into the thousand genomes data.')
    parser.add_argument('--data', type=str, dest="data",
                        help='Directory housing decompressed .vcf files 1-22,X')
    parser.add_argument('--populations', type=str, dest="pop",
                        help="Directory containing .txt lists of KGP sample IDs; the file name should correspond to the population.")
    parser.add_argument('--indexLevel', type=int, dest="level",
                        help="Every nth position will be indexed; a low number means a big index but faster lookup times. Probably a good balance is 100")
    
    args = parser.parse_args()
    
    print "Loading populations..."
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
    
    print "Loading positions..."
    # {chrN:{int(position/args.level):byte offset in file}}
    positions = {}
    # {chrN:file name}
    files = {}
    
    for dirname, dirnames, filenames in os.walk(args.data):
        for filename in filenames:
            if filename.endswith('.vcf') and 'chr' in filename:
                chrname = filename[filename.find("chr"):]
                chrname = chrname[:chrname.find(".")]
                fullpath = os.path.join(dirname,filename)
                files[chrname]=fullpath
                positions[chrname]={}
                
                infile = open(fullpath,'r')
                lastIndex = -1
                lastLine = 0
                for line in infile:
                    stripline = line.strip()
                    if len(stripline) <= 1 or line.startswith("##"):
                        lastLine = infile.tell()
                        continue
                    elif line.startswith("#"):
                        lastLine = infile.tell()
                        columns = line.strip().split('\t')
                        for i,c in enumerate(columns[9:]):
                            if reversePopulations.has_key(c):
                                if populations[reversePopulations[c]][c] == None:
                                    populations[reversePopulations[c]][c] = i+9
                                elif populations[reversePopulations[c]][c] != i+9:
                                    raise Exception("KGP files have %s in different columns: %i, %i" % (c,i+9,populations[reversePopulations[c]][c]))
                    else:
                        columns = stripline.split('/t')
                        currentPosition = int(columns[1])
                        currentIndex = currentPosition/args.level
                        if currentIndex > lastIndex:
                            positions[chrname][currentIndex]=lastLine
                        lastLine = infile.tell()
                infile.close()
    
    filestr = os.path.join(args.data,"KGP.index")
    outfile = open(filestr,'wb')
    pickle.dump(args.level,outfile)
    print "Writing populations..."
    pickle.dump(populations,outfile)
    print "Writing positions..."
    pickle.dump(positions,outfile)
    pickle.dump(files,outfile)
    outfile.close()
