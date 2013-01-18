#!/usr/bin/env python
import argparse, sys, os, gzip, math
from genome_utils import kgpInterface, countingDict, parsePopulations

class allStats:
    AF = 0
    Carriage = 1
    Samples_w_calls = 2
    
    STAT_NAMES=['AF',
                'Carriage',
                'SamplesWcalls']
    
    @staticmethod
    def calculate(stat,vcfLine,vcfIndices,kgpLine,kgpIndices,alleles):
        if kgpLine == None:
            kgpIndices = []
        if stat == allStats.AF:
            return allStats.calcAF(vcfLine,vcfIndices,kgpLine,kgpIndices,alleles)
        elif stat == allStats.Carriage:
            return allStats.calcCarriage(vcfLine,vcfIndices,kgpLine,kgpIndices,alleles)
        elif stat == allStats.Samples_w_calls:
            return allStats.calcSamples_w_calls(vcfLine,vcfIndices,kgpLine,kgpIndices,alleles)
        else:
            raise Exception("Unknown statistic: %s" % str(stat))
    
    @staticmethod
    def calcAF(vcfLine,vcfIndices,kgpLine,kgpIndices,alleles):
        if len(vcfIndices) > 0:
            vcfLine.extractAlleles()
            vcfLine.extractGenotypes(vcfIndices)
        if len(kgpIndices) > 0:
            kgpLine.extractAlleles()
            kgpLine.extractGenotypes(kgpIndices)
        
        count = 0.0
        matches = countingDict()
        for i in vcfIndices:
            allele0 = vcfLine.genotypes[i][0]
            allele1 = vcfLine.genotypes[i][1]
            if allele0 != None:
                count += 2.0
                matches[allele0] += 1
                matches[allele1] += 1
        for i in kgpIndices:
            allele0 = kgpLine.genotypes[i][0]
            allele1 = kgpLine.genotypes[i][1]
            if allele0 != None:
                count += 2.0
                matches[allele0] += 1
                matches[allele1] += 1
        
        if count == 0.0:
            return [float('Inf') for a in alleles]
        else:
            return [matches.get(i,0)/count for i,a in enumerate(alleles)]
        
    @staticmethod
    def calcCarriage(vcfLine,vcfIndices,kgpLine,kgpIndices,alleles):
        if len(vcfIndices) > 0:
            vcfLine.extractAlleles()
            vcfLine.extractGenotypes(vcfIndices)
        if len(kgpIndices) > 0:
            kgpLine.extractAlleles()
            kgpLine.extractGenotypes(kgpIndices)
        
        counts = countingDict()
        for i in vcfIndices:
            allele0 = vcfLine.genotypes[i][0]
            allele1 = vcfLine.genotypes[i][1]
            
            if allele0 != None:
                if allele0 == allele1:
                    counts[allele0] += 1
                else:
                    counts[allele0] += 1
                    counts[allele1] += 1
        for i in kgpIndices:
            allele0 = kgpLine.genotypes[i][0]
            allele1 = kgpLine.genotypes[i][1]
            
            if allele0 != None:
                if allele0 == allele1:
                    counts[allele0] += 1
                else:
                    counts[allele0] += 1
                    counts[allele1] += 1
        return [counts.get(i,0) for i,a in enumerate(alleles)]
        
    @staticmethod
    def calcSamples_w_calls(vcfLine,vcfIndices,kgpLine,kgpIndices,alleles):
        if len(vcfIndices) > 0:
            vcfLine.extractGenotypes(vcfIndices)
        if len(kgpIndices) > 0:
            kgpLine.extractGenotypes(kgpIndices)
        count = 0
        for i in vcfIndices:
            if vcfLine.genotypes[i][0] != None:
                count += 1
        for i in kgpIndices:
            if kgpLine.genotypes[i][0] != None:
                count += 1
        return count

count = 0
def tick():
    global count
    print "Calculating: {0}%\r".format(count),
    count += 1
    return True

def parseVcfHeader(path, outpath=None, popFile=""):
    # This is pretty ugly - it should probably be broken into several parts and called separately, but I didn't have time to do this right...
    # Essentially this pulls a lot of crap out of the vcf header, possibly in conjunction with a population file (see KGP_populations.txt or
    # a VCF Cleaner log file)
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
                    if i in columns[9:]:
                        populationIndices[p].append(columns.index(i)-9)
                    else:
                        populationIndices[p].append(i)   # this only happens if one of the names in the population file doesn't exist in the .vcf file... I assume it's someone from KGP and raise an Exception later if they're not
            infile.close()
            if outfile != None:
                return (outfile,infoTags,headerline,populations,populationIndices)
            return (infoTags,headerline,populations,populationIndices)
        else:
            raise Exception("Missing a header line or something else is wrong...")
    infile.close()

def run(args, tickFunction=tick, numTicks=100):
    kgp = kgpInterface(args.data,sys.path[0] + "/KGP_populations.txt")
    outfile,takenTags,headerline,myPopulations,myPopulationIndices = parseVcfHeader(args.infile,args.outfile,args.popFile)
    
    statsToCalculate = {}   # {INFO tag : (allStats.statistic,targetPop,backgroundPop,"ASC"/"DEC",REF/ALT hack: True/False,backTag))}
    alleleOrders = {}
    
    def storeCalcDetails(stat,calculation):
        if not len(calculation) > 0:
            raise Exception('Must specify a target population!')
        target = calculation[0]
        background = calculation[1] if len(calculation) > 1 else None
        direction = calculation[2] if len(calculation) > 2 else 'ASC'
        hack = len(calculation) > 3 and calculation[3].strip().lower().startswith('t')
        tag = "%s_%s_" % (target,allStats.STAT_NAMES[stat])
        if background == None:
            backTag = "ALT"
            tag += backTag
            infoLine = "##INFO=<ID=%s,Number=A,Type=Float,Description=\"calcStats.py: %s for the %s population\">\n" % (tag,allStats.STAT_NAMES[stat],target)
        else:
            backTag = "%s_%s_AO" % (direction,background)
            temp = backTag
            dupNumber = 2
            while backTag in takenTags or backTag in myPopulations.iterkeys() or backTag in kgp.populations.iterkeys():
                backTag = temp + str(dupNumber)
                dupNumber += 1
            alleleOrders[backTag] = (direction,background)
            tag += backTag
            if hack:
                tag += "_rHack"
            infoLine = "##INFO=<ID=%s,Number=.,Type=Float,Description=\"calcStats.py: %s for the %s population, with alleles ordered by %s AF in the %s population (%s).%s\">\n" % (tag,
                       allStats.STAT_NAMES[stat],
                       target,
                       "ascending" if direction == 'ASC' else "descending",
                       background,
                       backTag,
                       " When %s has no data, the REF/ALT allele order is used." % background if hack else "")
        dupNumber = 2
        temp = tag
        while tag in takenTags:
            tag = temp + str(dupNumber)
            dupNumber += 1
        statsToCalculate[tag] = (stat,target,background,direction,hack,backTag)
        return infoLine
    
    if args.calculate_AF != None:
        for calculation in args.calculate_AF:
            outfile.write(storeCalcDetails(allStats.AF,calculation))
    if args.calculate_Carriage != None:
        for calculation in args.calculate_Carriage:
            outfile.write(storeCalcDetails(allStats.Carriage,calculation))
    if args.calculate_Samples_w_calls != None:
        for calculation in args.calculate_Samples_w_calls:
            assert len(calculation) == 1
            outfile.write(storeCalcDetails(allStats.Samples_w_calls,calculation))
    
    for popTag,(direction,background) in alleleOrders.iteritems():
        outfile.write("##INFO=<ID=%s,Number=.,Type=String,Description=\"calcStats.py: All observed alleles for each locus, ordered by %s AF in the %s population.\">\n" % (popTag,
                      "ascending" if direction == 'ASC' else "descending",
                      background))
    outfile.write(headerline)
    
    def getPopIndices(pop):
        if myPopulationIndices.has_key(pop):
            vcfIndices = []
            kgpIndices = []
            for i in myPopulationIndices[pop]:
                if isinstance(i,str):
                    if not i in kgp.header[9:]:
                        raise Exception("Unknown sample (not in your .vcf or the KGP): %s" % i)
                    i = kgp.header.index(i)-9
                    kgpIndices.append(i)
                else:
                    vcfIndices.append(i)
        else:
            vcfIndices = []
            kgpIndices = kgp.populationIndices[pop]
        return (vcfIndices,kgpIndices)
    
    for vcfLine,kgpLine in kgp.iterateVcf(args.infile,tickFunction=tickFunction,numTicks=numTicks):
        # first get the allele orders we need, add them as INFO fields
        alleleLists = {}    # popTag : []
        vcfLine.extractAlleles()
        vcfLine.extractInfo()
        if kgpLine != None:
            kgpLine.extractAlleles()
        
        for popTag,(direction,background) in alleleOrders.iteritems():
            vcfIndices,kgpIndices = getPopIndices(background)
            tempAlleles = set(vcfLine.alleles)
            if kgpLine != None:
                tempAlleles.update(kgpLine.alleles)
            tempAlleles = list(tempAlleles)
            tempFreqs = allStats.calculate(allStats.AF,vcfLine,vcfIndices,kgpLine,kgpIndices,tempAlleles)
            if len(tempFreqs) < 1 or math.isinf(tempFreqs[0]):
                vcfLine.info[popTag] = "."
                alleleLists[popTag] = None
            else:
                if direction == 'ASC':
                    alleleLists[popTag] = sorted(tempAlleles,key=lambda i:tempFreqs[tempAlleles.index(i)])
                else:
                    alleleLists[popTag] = sorted(tempAlleles,key=lambda i:tempFreqs[tempAlleles.index(i)],reverse=True)
                vcfLine.info[popTag] = ",".join(alleleLists[popTag])
        
        # now calculate based on those allele orders
        for tag,(stat,target,background,direction,hack,backTag) in statsToCalculate.iteritems():
            vcfIndices,kgpIndices = getPopIndices(target)
            if backTag == 'ALT':
                alleles = vcfLine.alleles
            else:
                alleles = alleleLists[backTag]
            if alleles == None:
                if hack:
                    alleles = vcfLine.alleles
                else:
                    vcfLine.info[tag] = "."
                    continue
            result = allStats.calculate(stat,vcfLine,vcfIndices,kgpLine,kgpIndices,alleles)
            if isinstance(result,list):
                result = ",".join([str(r) for r in result])
            else:
                result = str(result)
            vcfLine.info[tag] = result
        
        outfile.write(str(vcfLine))
    
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add some calculated statistics to a .vcf file\'s INFO column. For each --calculate parameter, '+
                        'supply the population name (the header line) from either your samples (--populations) or the 1000 Genomes Project (KGP_populations.txt)')
    parser.add_argument('--in', type=str, dest="infile", required=True,
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile", required=True,
                        help='output .vcf file')
    parser.add_argument('--data', type=str, dest="data", required=True,
                        help='Path to directory containing 1000 genomes .vcf.gz files 1-22,X,Y.')
    parser.add_argument('--populations', type=str, dest="popFile", nargs="?", const="", default="",
                        help='Population file describing samples in your .vcf file. See the VCF Cleaner help page (http://sci.utah.edu/~abigelow/vcfCleanerHelp.php#Populations) for details.')
    parser.add_argument('--calculate_AF', type=str, dest="calculate_AF", nargs="+", action="append",
                        help='Recalculates allele frequencies, may be used multiple times. At least one argument is required. The first argument should be the population in which to calculate allele frequencies. The second is '+
                        'the background population; per-allele values will be ordered by AF in this population. If omitted, the ALT allele order from --in is used. The third argument is "ASC" or "DEC", indicating '+
                        'whether the alleles are ordered in ascending or descending AF: if omitted, "ASC" is assumed. The fourth argument is "True" or "False"; if true, the REF/ALT order will be reused when '+
                        'a variant has no data in the background population. As this is technically a violation of nomenclature, if omitted "False" is assumed.')
    parser.add_argument('--calculate_Carriage', type=str, dest="calculate_Carriage", nargs="+", action="append",
                        help='Counts the number of individuals who have at least one copy of an allele, may be used multiple times. At least one argument is required. The first argument should be the population in which to calculate carriage of each allele. The second is '+
                        'the background population; per-allele values will be ordered by AF in this population. If omitted, the ALT allele order from --in is used. The third argument is "ASC" or "DEC", indicating '+
                        'whether the alleles are ordered in ascending or descending AF: if omitted, "ASC" is assumed. The fourth argument is "True" or "False"; if true, the REF/ALT order will be reused when '+
                        'a variant has no data in the background population. As this is technically a violation of nomenclature, if omitted "False" is assumed.')
    parser.add_argument('--calculate_Samples_w_calls', type=str, dest="calculate_Samples_w_calls", nargs="+", action="append",
                        help='Calculates the number of samples with calls, may be used multiple times. Exactly one argument is required. The argument should be the population in which to count samples with calls.')
    
    args = parser.parse_args()
    run(args)