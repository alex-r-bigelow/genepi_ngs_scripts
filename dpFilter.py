#!/usr/bin/env python
import argparse, math

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Applies "DPfilter" to variants that are "5 or 6 sigma from the '+
                                     'mean coverage across all samples" (see GATK best practice v4')
    parser.add_argument('--stddev', type=int, dest="stddev",
                        help='number of standard deviations below which a variant will PASS.')
    parser.add_argument('--in', type=str, dest="infile",
                        help='input .vcf file')
    parser.add_argument('--out', type=str, dest="outfile",
                        help='output .vcf file')
    
    args = parser.parse_args()
    
    sum = 0
    squaresum = 0
    count = 0
    print "Calculating DP standard deviation"
    print "Reading..."
    infile = open(args.infile,'r')
    for line in infile:
        if len(line.strip()) <= 1 or line.startswith("#"):
            continue
        columns = line.split()
        info = columns[7].split(";")
        for i in info:
            if i.startswith('DP'):
                dp = i.split("=")
                depth = int(dp[1])
                count += 1
                sum += depth
                squaresum += depth**2
                break
    infile.close()
    
    mean = float(sum)/count
    sigma = math.sqrt(float(squaresum)/count)
    print "Mean: %f Sigma: %f" % (mean,sigma)
    
    print "Writing..."
    wrotePragma = False
    
    outfile = open(args.outfile,'w')
    infile = open(args.infile,'r')
    for line in infile:
        if len(line.strip()) <= 1:
            continue
        if line.startswith("#"):
            if not wrotePragma and line.startswith("##FILTER"):
                # Stick our pragma line in before the other filters
                outfile.write("##FILTER=<ID=DP Filter,Description=\"Depth Filter; DP >= %i sigma\">\n")
                wrotePragma = True
            outfile.write(line)
        else:
            columns = line.split()
            filters = columns[6].split(";")
            info = columns[7].split(";")
            for i in info:
                if i.startswith('DP'):
                    dp = i.split("=")
                    depth = int(dp[1])
                    stddevs = abs(depth-mean)/sigma
                    if stddevs >= args.stddev:
                        if 'PASS' in filters:
                            filters.remove('PASS')
                        filters.append('DP Filter')
                    break
            columns[6] = ";".join(filters)
            outfile.write('\t'.join(columns)+"\n")
    infile.close()
    outfile.close()
    print "Done"