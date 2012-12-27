#!/bin/bash

# downloads the compressed KGP .vcf files to the current directory

for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	curl -O http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr$c.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
done