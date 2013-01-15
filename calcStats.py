#!/usr/bin/env python
import argparse, sys, os, gzip
from genome_utils import kgpInterface, countingDict, parsePopulations

class allStats:
    AF = 0
    Carriage = 1
    Samples_w_calls = 2
    
    STAT_NAMES=['AF',
                'Carriage',
                'Samples w/Calls']
    
    @staticmethod
    def calculate(stat,line,indices):
        if stat == allStats.AF:
            return allStats.calcAF(line, indices)
        elif stat == allStats.Carriage:
            return allStats.calcCarriage(line, indices)
        elif stat == allStats.Samples_w_calls:
            return allStats.calcSamples_w_calls(line, indices)
        else:
            raise Exception("Unknown statistic: %s" % str(stat))
    
    @staticmethod
    def calcAF(line, indices):
        line.extractAlleles()
        line.extractGenotypes(indices)
        
        count = 0.0
        matches = countingDict()
        for i in indices:
            allele0 = line.genotypes[i][0]
            allele1 = line.genotypes[i][1]
            if allele0 != None:
                count += 2.0
                matches[allele0] += 1
                matches[allele1] += 1
        
        if count == 0.0:
            return [float('Inf') for a in line.alleles]
        else:
            return [matches.get(i,0)/count for i,a in enumerate(line.alleles)]
        
    @staticmethod
    def calcCarriage(line, indices):
        line.extractAlleles()
        line.extractGenotypes(indices)
        
        counts = countingDict()
        for i in indices:
            allele0 = line.genotypes[i][0]
            allele1 = line.genotypes[i][1]
            
            if allele0 != None:
                if allele0 == allele1:
                    counts[allele0] += 1
                else:
                    counts[allele0] += 1
                    counts[allele1] += 1
        return [counts.get(i,0) for i,a in enumerate(line.alleles)]
        
    @staticmethod
    def calcSamples_w_calls(line, indices):
        line.extractGenotypes(indices)
        count = 0
        for i in indices:
            if line.genotypes[i][0] != None:
                count += 1
        return count

def tick():
    print ".",

def parseVcfHeader(path, outpath=None, popFile=""):
    if popFile != "":
        populations = parsePopulations(args.popFile)[0]
    else:
        popName = os.path.split(path)[1]
        populations = {popName:[]}
    populationIndices = {}
    
    infoTags = set()
    
    if path.endswith('.gz'):
        infile = gzip.open(path,'rb')
    else:
        infile = open(path,'rb')
    
    outfile = None
    if outpath != None:
        if outpath.endswith('.gz'):
            outfile = gzip.open(outpath,'wb')
        else:
            outfile = open(outpath, 'wb')
    
    for line in infile:
        if len(line) <= 1:
            continue
        elif line.startswith("##"):
            if line.startswith("##INFO"):
                newTag = line[line.find("ID=")+3:]
                newTag = newTag[:newTag.find(',')]
                infoTags.add(newTag)
            if outpath != None:
                outfile.write(line)
        elif line.startswith("#"):
            headerline = line
            columns = line.strip().split('\t')
            if popFile == "":
                populations[popName] = columns[9:]
            for p,individuals in populations.iteritems():
                populationIndices[p] = []
                for i in individuals:
                    populationIndices[p].append(columns.index(i)-9)
            infile.close()
            if outfile != None:
                return (outfile,infoTags,headerline,populations,populationIndices)
            return (infoTags,headerline,populations,populationIndices)
        else:
            raise Exception("Missing a header line or something else is wrong...")
    infile.close()

def extractCalcDetails(tag,):

def run(args, tickFunction=tick):
    outfile,takenTags,headerline,myPopulations,myPopulationIndices = parseVcfHeader(args.infile,args.outfile,args.popFile)
    
    statsToCalculate = {}   # {INFO tag : (allStats.statistic,targetPop,backgroundPop,"ASC"/"DEC",REF/ALT hack: True/False))}
    alleleOrders = {}       # {backgroundPop : "ASC","DEC", or "BOTH"}
    
    def extractCalcDetails(stat,target,background,direction,hack):
        
    
    if args.calculate_AF != None:
        for calculation in args.calculate_AF:
            if not len(calculation) > 0:
                raise Exception('Must specify a target population!')
            else:
                if not len(calculation) > 1:
                    statsToCalculate[calculation[0]]
            
            
            dupNumber = 2
            tag = p + "_AF"
            while tag in takenTags:
                tag = p + "_AF" + str(dupNumber)
                dupNumber += 1
            statsToCalculate[tag] = (allStats.AF,p)
            outfile.write("##INFO=<ID=%s,Number=.,Type=Float,Description=\"calcStats.py: All allele frequencies in the %s population\">\n" % (tag,p))
    if args.calculate_MAF != None:
        for p in args.calculate_MAF:
            dupNumber = 2
            mafTag = "_MAF"
            if args.reorder_alleles == "NO_REORDERING":
                mafTag = "_pseudoMAF"
            tag = p + mafTag
            while tag in takenTags:
                tag = p + mafTag + str(dupNumber)
                dupNumber += 1
            statsToCalculate[tag] = (allStats.MAF,p)
            if args.reorder_alleles == "NO_REORDERING":
                outfile.write("##INFO=<ID=%s,Number=A,Type=Float,Description=\"calcStats.py: Minor allele frequency(ies) in the %s population. WARNING: these rely on the REF/ALT configuration of the calling pipeline; THEY ARE NOT TRUE MINOR ALLELE FREQUENCIES!\">\n" % (tag,p))
            else:
                outfile.write("##INFO=<ID=%s,Number=A,Type=Float,Description=\"calcStats.py: Minor allele frequency(ies) in the %s population\">\n" % (tag,p))
    if args.calculate_Sharing != None:
        for p in args.calculate_Sharing:
            dupNumber = 2
            tag = p + "_Sharing"
            while tag in takenTags:
                tag = p + "_Sharing" + str(dupNumber)
                dupNumber += 1
            statsToCalculate[tag] = (allStats.Sharing,p)
            outfile.write("##INFO=<ID=%s,Number=.,Type=Integer,Description=\"calcStats.py: Maximum sharing for every allele in the %s population\">\n" % (tag,p))
    if args.calculate_Minor_Sharing != None:
        for p in args.calculate_Minor_Sharing:
            dupNumber = 2
            mafTag = "_Minor_Sharing"
            if args.reorder_alleles == "NO_REORDERING":
                mafTag = "_pseudoMinor_Sharing"
            tag = p + mafTag
            while tag in takenTags:
                tag = p + mafTag + str(dupNumber)
                dupNumber += 1
            statsToCalculate[tag] = (allStats.Minor_Sharing,p)
            if args.reorder_alleles == "NO_REORDERING":
                outfile.write("##INFO=<ID=%s,Number=A,Type=Float,Description=\"calcStats.py: Maximum sharing for every minor allele in the %s population. WARNING: these rely on the REF/ALT configuration of the calling pipeline; THEY ARE NOT TRUE MINOR ALLELE FREQUENCIES!\">\n" % (tag,p))
            else:
                outfile.write("##INFO=<ID=%s,Number=A,Type=Float,Description=\"calcStats.py: Maximum sharing for every minor allele in the %s population\">\n" % (tag,p))
    if args.calculate_Samples_w_calls != None:
        for p in args.calculate_Samples_w_calls:
            dupNumber = 2
            tag = p + "_Samples_w_calls"
            while tag in takenTags:
                tag = p + "_Samples_w_calls" + str(dupNumber)
                dupNumber += 1
            statsToCalculate[tag] = (allStats.Samples_w_calls,p)
            outfile.write("##INFO=<ID=%s,Number=1,Type=Integer,Description=\"calcStats.py: Number of samples with calls in the %s population\">\n" % (tag,p))
    
    outfile.write(headerline)
    
    kgp = kgpInterface(args.data,sys.path[0] + "/KGP_populations.txt")
    
    for line,kgpLine in kgp.iterateVcf(args.infile,tickFunction=tickFunction,numTicks=1000):
        if args.reorder_alleles != "NO_REORDERING":
            if myPopulationIndices.has_key(args.reorder_alleles):
                temp = allStats.calcAF(line, myPopulationIndices[args.reorder_alleles])
                # replace any undefined frequencies with zero, as this population never sees any alleles!
                if float('Inf') in temp:
                    targetFreqs = {}
                else:
                    targetFreqs = dict(zip(line.alleles,temp))
            elif kgp.populationIndices.has_key(args.reorder_alleles):
                if kgpLine == None:
                    targetFreqs = {}
                else:
                    temp = allStats.calcAF(kgpLine, kgp.populationIndices[args.reorder_alleles])
                    # replace any undefined frequencies with zero, as this population never sees any alleles!
                    if float('Inf') in temp:
                        targetFreqs = {}
                    else:
                        targetFreqs = dict(zip(kgpLine.alleles,temp))
            else:
                raise Exception("Unknown population: %s" % args.reorder_alleles)
            line.reorderAlleles(targetFreqs)
            if kgpLine != None:
                kgpLine.reorderAlleles(targetFreqs)
        
        line.extractInfo()
        for tag,(stat,pop) in statsToCalculate.iteritems():
            if myPopulationIndices.has_key(pop):
                result = allStats.calculate(stat,line,myPopulationIndices[pop])
            elif kgp.populationIndices.has_key(pop):
                if kgpLine == None:
                    result = "."
                else:
                    result = allStats.calculate(stat,kgpLine,kgp.populationIndices[pop])
            else:
                raise Exception("Unknown population: %s" % pop)
            
            if isinstance(result,list):
                result = ",".join([str(r) for r in result])
            else:
                result = str(result)
            line.info[tag] = result
        
        outfile.write(str(line))
    
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add some calculated statistics to a .vcf file\'s INFO column. For each --calculate parameter, '+
                        'supply the population name (the header line) from either your samples (--populations) or the 1000 Genomes Project (KGP_populations.txt)')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    parser.add_argument('--data', type=str, dest="data",
                        help='Path to directory containing 1000 genomes .vcf.gz files 1-22,X,Y.')
    parser.add_argument('--populations', type=str, dest="popFile", nargs="?", const="", default="",
                        help='Tab-delimited .txt file containing populations in your .vcf file. A header is required, (each header is the population name), and each row under the header is a sample ID from the .vcf file indicating that that sample is a member of that population. If not supplied, all samples will be used and the population name will be the file name (e.g. "myFile.vcf"). '+
                        'See KGP_populations.txt for an example of how to format this file. Be careful to not reuse any of the column headers in KGP_populations.txt for your own data.')
    parser.add_argument('--calculate_AF', type=str, dest="calculate_AF", nargs="+", action="append",
                        help='Recalculates allele frequencies, may be used multiple times. First argument should be the population in which to calculate allele frequencies. The second is '+
                        'the background population: if omitted, the REF/ALT configuration of --in is used. The third argument is "ASC" or "DEC", indicating the '+
                        'order of the alleles: if omitted, "ASC" is assumed. The fourth argument is "True" or "False"; if true, the REF/ALT order will be reused when '+
                        'a variant has no data in the background population. As this is technically a violation of nomenclature, if omitted "False" is assumed.')
    parser.add_argument('--calculate_Carriage', type=str, dest="calculate_Carriage", nargs="+", action="append",
                        help='Counts the number of individuals who have at least one copy of an allele, may be used multiple times. First argument should be the population in which to calculate allele frequencies. The second is '+
                        'the background population: if omitted, the REF/ALT configuration of --in is used. The third argument is "ASC" or "DEC", indicating the '+
                        'order of the alleles: if omitted, "ASC" is assumed. The fourth argument is "True" or "False"; if true, the REF/ALT order will be reused when '+
                        'a variant has no data in the background population. As this is technically a violation of nomenclature, if omitted "False" is assumed.')
    parser.add_argument('--calculate_Samples_w_calls', type=str, dest="calculate_Samples_w_calls", nargs="+", action="append",
                        help='Calculates the number of samples with calls, may be used multiple times. First argument should be the population in which to calculate allele frequencies. The second is '+
                        'the background population: if omitted, the REF/ALT configuration of --in is used. The third argument is "ASC" or "DEC", indicating the '+
                        'order of the alleles: if omitted, "ASC" is assumed. The fourth argument is "True" or "False"; if true, the REF/ALT order will be reused when '+
                        'a variant has no data in the background population. As this is technically a violation of nomenclature, if omitted "False" is assumed.')
    
    args = parser.parse_args()
    run(args)