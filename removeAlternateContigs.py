#!/usr/bin/env python
import argparse

if __name__ == '__main__':
    chromosomesToKeep = set(['chr1',
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
                             'chrY'])
    
    parser = argparse.ArgumentParser(description='Leaves only chromosomes 1-22,X,Y in a .vcf file')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    
    args = parser.parse_args()
    
    infile = open(args.infile,'r')
    outfile = open(args.outfile,'w')
    for line in infile:
        if line.startswith("#"):
            if 'contig' in line:
                id = line[line.find('ID=')+3:]
                before = line[:line.find('ID=')+3]
                after = id[id.find(','):]
                id = id[:id.find(',')]
                
                if not id.startswith('chr'):
                    id = 'chr' + id
                if id in chromosomesToKeep:
                    outfile.write(before + id + after)
            else:
                outfile.write(line)
                continue
        columns = line.strip()
        if len(columns) <= 1:
            continue
        columns = columns.split('\t')
        chr = columns[0]
        if not chr.startswith('chr'):
            chr = "chr" + chr
        if chr in chromosomesToKeep:
            outfile.write("%s\t%s\n"%(chr,'\t'.join(columns[1:])))
    infile.close()
    outfile.close()