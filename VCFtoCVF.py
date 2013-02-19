#!/usr/bin/env python
import argparse, datetime, math, sys
from genome_utils import standardizeChromosome, vcfLine, countingDict, infoDetails, MAX_INFO_STRINGS

def run(args):
    
    separateInfoFields = args.separate_info_fields.strip().lower() == "true"
    countSeparate = args.count_separate.strip().lower() == "true"
    
    ignoreFields = args.ignore_fields
    if ignoreFields == None:
        ignoreFields = []
    ignoreFields = set(ignoreFields)
    
    posLengthWarned = False
    allChrs = []
    positions = []
    alleleColumn = infoDetails("Ref/Alt", 1, False)
    qualColumn = infoDetails("QUAL", 1, False)
    filterColumn = infoDetails("FILTER", args.max_strings, False)
    infoFields = {"Ref/Alt":alleleColumn,"QUAL":qualColumn,"FILTER":filterColumn}
    # TODO: get the numeric ranges, all valid categorical values
    
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
                if infoFields.has_key(newTag):
                    raise Exception("Duplicate INFO ID or use of reserved ID:\t%s" % newTag)
                infoFields[newTag] = infoDetails(newTag, args.max_strings, countSeparate)
                if newTag in ignoreFields:
                    infoFields[newTag].maxedOut = True
            elif temp.startswith("##filter"):
                newTag = line[temp.find("id=")+3:]
                newTag = newTag[:newTag.find(',')]
                filterColumn.addCategory(newTag, 0)
            elif temp.startswith("##contig"):
                chrom = line[temp.find("id=")+3:]
                chrom = newTag[:chrom.find(',')]
                chrom = standardizeChromosome(chrom)
                chrLength = line[temp.find("length=")+3:]
                chrLength = chrLength[:chrLength.find(',')]
                
                allChrs.append(chrom)
                positions.append((0,int(chrLength)))
            else:
                # a sneaky way of freezing the filter column; if other filters are added (without a .vcf pragma line) or we aren't separating info fields,
                # other strings will make this column max out early
                filterColumn.maxCategories = len(filterColumn.categories[0])
        else:
            line = vcfLine(line.split('\t'))
            line.extractChrAndPos()
            
            if not line.chromosome not in allChrs:
                allChrs.append(line.chromosome)
                positions.append((0,0))
            chrIndex = allChrs.index(line.chromosome)
            if line.position > positions[chrIndex][1]:
                positions[chrIndex] = (0,line.position)
                if not posLengthWarned:
                    sys.stderr.write('WARNING: Either ##contig pragma lines aren\'t supplied in your .vcf file or a variant has a position beyond the length of a chromosome.')
                    sys.stderr.write(' In either case, be aware that chromosome lengths in the .cvf file may not be accurate.')
                    posLengthWarned = True
            
            line.extractAlleles()
            alleles = line.alleles
            if not separateInfoFields:
                alleles = ",".join(alleles)
            alleleColumn.addArbitraryValue(line.alleles)
            
            line.extractQual()
            qualColumn.addArbitraryValue(line.qual)
            
            line.extractFilters()
            filters = line.filters
            if not separateInfoFields:
                filters = ",".join(filters)
            filterColumn.addArbitraryValue(filters)
            
            line.extractInfo()
            for k,v in line.info.iteritems():
                if not infoFields.has_key(k):
                    raise Exception("Missing ##INFO pragma for: %s" % k)
                if separateInfoFields:
                    v = ",".split(v)
                infoFields[k].addArbitraryValue(v)
    infile.close()
    
    print "Creating file..."
    outfile = open(args.outfile, 'w')
    
    outfile.write("##\t%s created from %s on %s\n" % (args.outfile,args.infile,str(datetime.datetime.now())))
    outfile.write("#\tChromosome\tCHR\t%s\n" % ("\t".join(allChrs)))
    outfile.write("#\tPosition\tPOS\t%s\n" % ("\t".join(["(%i,%i)" % p for p in positions])))
    outfile.write("#\tID\tID\n")
    
    headers = []
    fieldOrder = sorted(infoFields.iterkeys())
    for f in fieldOrder:
        pragmas = infoFields[f].getPragmas()
        for p in pragmas:
            outfile.write(p + "\n")
            h = p.split("\t")[1]
            headers.append(h)
    
    outfile.write('Chromosome\tPosition\tID\t%s\n' % ("\t".join(headers)))
    
    infile = open(args.infile,'r')
    for line in infile:
        line = line.strip()
        if len(line) <= 1 or line.startswith("#"):
            continue
        line = vcfLine(line.split('\t'))
        line.extractChrAndPos()
        line.extractInfo()
        line.extractAlleles()
        line.info["Ref/Alt"] = line.alleles
        line.extractQual()
        line.info["QUAL"] = str(line.qual)
        line.extractFilters()
        line.info["FILTER"] = line.filters
        outfile.write("%s\t%i\t%s" % (line.chromosome,line.position,line.id))
        for f in fieldOrder:
            values = line.info[f]
            if isinstance(values,list):
                if separateInfoFields:
                    values = "\t".join(values)
                else:
                    values = ",".join(values)
            outfile.write("\t%s" % values)
        outfile.write("\n")
    infile.close()
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates a .cvf file from the CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO fields of a .vcf file')
    parser.add_argument('--in', type=str, dest="infile", required=True,
                        help='Path to .vcf file')
    parser.add_argument('--out', type=str, dest="outfiles", required=True,
                        help='Path to .cvf file')
    parser.add_argument('--max_strings', type=int, dest="max_strings", nargs="?", const=MAX_INFO_STRINGS, default=MAX_INFO_STRINGS,
                        help='Maximum number of strings a categorical INFO field can have before it\'s automatically marked as IGNORE. If zero or negative, no limit is enforced. Default is %i.' % MAX_INFO_STRINGS)
    parser.add_argument('--separate_info_fields', type=str, dest="separate_info_fields", nargs="?", const="True", default="False",
                        help='When an INFO field has multiple comma-delimited values, split them into separate columns.')
    parser.add_argument('--count_separate', type=str, dest="count_separate", nargs="?", const="True", default="False",
                        help='If --separate_info_fields is supplied, this option will count each column\'s ranges and categorical values separately.')
    parser.add_argument('--ignore', type=str, dest="ignore_fields", nargs="+",
                        help='Flag specific columns as IGNORE; this overrides --max_strings.')
    
    args = parser.parse_args()
    run(args)