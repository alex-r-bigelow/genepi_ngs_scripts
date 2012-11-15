genepi_ngs_scripts
==================

Scripts, tools for working with next-generation sequencing data

best_practice_v4_exomes.sh, best_practice_v4_genomes.sh
-------------------------------------------------------
BWA/GATK calling pipeline per GATK Best Practice v4

runAnnotationTools.sh
---------------------
Runs my dpFilter.py, calcStats.py programs, as well as SnpEff, VAAST, and ANNOVAR

dpFilter.py
-----------
Per the Best Practice hard filters, a DP filter should be applied if the DP exceeds 5 or 6 sigma

calcStats.py
------------
An example script that calculates the following for every variant:
 - sharing for every vcf allele
 - allele frequency for every vcf allele
 - thousand genomes' allele frequency for every vcf allele
Also prints out how many variants would be left if certain filters were applied

This script is prone to rapid change

index_kgp.py
------------
Library for querying thousand genotypes' .vcf files