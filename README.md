genepi_ngs_scripts
==================

Scripts, tools for working with next-generation sequencing data. For a very simple example of how to use my library in your own scripts, see sampleParse.py

best_practice_v4.sh, example.sh
-------------------------------
BWA/GATK calling pipeline per GATK Best Practice v4 (http://www.broadinstitute.org/gatk/guide/topic?name=best-practices). Options aren't yet comprehensive enough to cover the whole spec: for now, this assumes the "Best" post-processing strategy and uses hard filters (I only work with a few samples at a time).

To use this, you should make a copy of example.sh (named appropriately) for each run that you do, tweaking its parameters as needed.

This also makes use of:

- dpFilter.py:
  Per the Best Practice hard filters, a DP filter should be applied if the DP exceeds 5 or 6 sigma, but as far as I could tell, there's no utility to do this in GATK

- buildVAASTreference.py:
  This is a script that attempts to tweak the GATK bundle's .fasta reference genome to be compatible with VAAST (technically they give you the same reference genome, but I'm OCD :) ). As VAAST is still under heavy development, this works only about half of the time.

- removeAlternateContigs.py:
  For VAAST compatibility, strips out variants in a .vcf file that are not in chromosomes 1-22,X, or Y

vcfCleaner.py
-------------
A GUI front end to scripts that can manipulate/clean the results of the pipeline. The GUI is not quite ready, but each script can run independently:

- sort.py:
  Should be run on any file that is to be fed to any of these other scripts (sorts chromosomes, positions in 1-22,X,Y. All other chromosomes (chrUn, MT, etc) come after in alphabetic order)

- addBEDtoVCF.py:
  Adds per-feature scores in a .bed file to every intersecting variant in a .vcf file

- addCSVtoVCF.py:
  Adds per-variant scores in a .csv file to every variant in a .vcf file; supports three modes for matching rows: exact match, nearest neighbor, and interpolation

- calcStats.py:
  Really does two things: creates additional/alternate allele orderings (see http://sci.utah.edu/~abigelow/vcfCleanerHelp.php#Calculating for an explanation), and calculates additional statistics for those allele orders

- cleanVCF.ph:
  Removes INFO fields from a .vcf file with an excessive number of categorical values

- filterVCF.py:
  Filters rows of a .vcf file per Python-syntax expressions

- VCFtoCVF.py:
  Converts a .vcf file to a .cvf file