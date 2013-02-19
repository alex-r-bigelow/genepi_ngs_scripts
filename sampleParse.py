from genome_utils import vcfLine

'''
A very basic sample script showing how to use vcfLine; takes a .vcf file as "input.vcf" and writes an "output.csv" file with each
genotype from the .vcf file.
'''

with open('input.vcf','rb') as infile, open('output.csv','wb') as outfile:
    individualIDs = None
    for line in infile:
        if len(line) <= 1:
            # Skip blank lines
            continue
        elif line.startswith('#'):
            # Skip header pragma lines
            if line.lower().startswith('#chrom'):
                # Get all the individual IDs from the header, write it to the file
                individualIDs = line.strip().split('\t')[9:]
                outfile.write(','.join(individualIDs))
                outfile.write('\n')
        else:
            assert individualIDs != None    # In a well-formed .vcf file, we'll have run across the header line before any data
            
            # Here is where my library comes in - we first build a vcfLine object with the columns in the .vcf file
            line = vcfLine(line.strip().split('\t'))
            
            # Ideally, these two steps would be performed automatically, but sometimes we might want
            # to skip them to save time in practice:
            
            # We have to first extract the alleles before I can reference line.alleles
            line.extractAlleles()
            # We have to first extract all genotypes before I can reference line.genotypes
            line.extractGenotypes()
            
            for i in xrange(len(individualIDs)):
                allele1,allele2,phased,attributes = line.genotypes[i]
                
                outfile.write(line.alleles[allele1])   # write the first allele letters - to write the number, just write allele1
                outfile.write("|" if phased else "/")
                outfile.write(line.alleles[allele2])   # write the second allele letters - to write the number, just write allele2
                if i == len(individualIDs) - 1:
                    outfile.write('\n') # we're on the last ID; write a newline
                else:
                    outfile.write('\t') # we're not on the last ID; write a tab
    infile.close()
    outfile.close()