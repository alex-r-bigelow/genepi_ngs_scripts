#!/bin/bash

BWA_DIR=/raid1/sequencing/apps/alignment/bwa-0.6.2
GATK_DIR=/raid1/sequencing/apps/post_processing/GenomeAnalysisTK-2.1-11-g13c0244
REF_DIR=/raid1/sequencing/reference/gatk_bundle/gatk_bundle_04_Oct_2012
DATA_DIR=/raid1/alex/sequencing/cll/data/genomes

NUM_CORES=8
RUN_JAVA="/usr/local/java/jdk1.6.0_22_x64/bin/java -Xmx16g"
