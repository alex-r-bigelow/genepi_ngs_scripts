#!/bin/bash
# To be run after one of the best_practice_v4 pipelines; will add an additional directory containing various
# annotations. The script will run dpFilter.py, calcStats.py, SnpEff, VAAST, and ANNOVAR

# ******** Helper Function ********
# this will make the whole script fail if any piece does
# (also waits for multiple parallel jobs to finish before
# continuing)
waitForJobs ()
{
	FAIL=0
	for job in `jobs -p`
	do
		wait $job || let "FAIL+=1"
	done
	if [[ $FAIL != 0 ]]
	then
		echo "A job failed!"
		echo "A job failed!" >&2
		exit 1
	fi
}

# ******** Setup ********
# <<<< edit these values to point to the appropriate programs, reference files >>>>
GATK_DIR=/raid1/sequencing/apps/post_processing/GenomeAnalysisTK-2.1-11-g13c0244
VAAST_DIR=/raid1/sequencing/apps/analysis/VAAST/VAAST_trunk/bin
SNPEFF_DIR=/raid1/sequencing/apps/annotation/snpEff_3_1
ANNOVAR_DIR=/raid1/sequencing/apps/annotation/annovar

KGP_DIR=/raid1/sequencing/reference/background/KGP/
KGP_DATA_DIR=$KGP_DIR/compressed_vcfs
KGP_POP_DIR=$KGP_DIR/populationLists
REF_DIR=/raid1/sequencing/reference/gatk_bundle/gatk_bundle_04_Oct_2012

REF_FASTA=$REF_DIR/human_g1k_v37.fasta

# <<<< adjust these two values to match the number of cores, GB of ram on the machine >>>>
NUM_CORES=12
MAX_MEM=64

MAX_BYTE_MEM=$(( $MAX_MEM*1024*1024*1024 ))
MAX_SAMPLE_MEM=$(( $MAX_BYTE_MEM/$NUM_SAMPLES ))
RUN_JAVA="/usr/local/java/jdk1.6.0_22_x64/bin/java -Xmx$MAX_SAMPLE_MEM"

$RUN_JAVA >/dev/null 2>/dev/null &
waitForJobs
# silently test if we have enough memory... but fail early if we don't

# <<<< these only need to be run once >>>>
python index_kgp.py --data $KGP_DATA_DIR --populations $KGP_POP_DIR --indexLevel 64
#$RUN_JAVA -jar $SNPEFF_DIR/snpEff.jar download -v hg19 &
waitForJobs

# <<<< reset this for a different project >>>>
TARGET_DIR=/raid1/alex/sequencing/cll/runs/exome_22_Oct_2012
TARGET_VCFS=(`ls $TARGET_DIR/calls/all.*.filtered.vcf`)

# I use this if statement to comment stuff out
#if [[ 0 == 1 ]]
#then
# ******** Actual jobs start here ********
rm -rf $TARGET_DIR/annotation
mkdir $TARGET_DIR/annotation

for s in ${TARGET_VCFS[*]}
do
	echo $s
done
