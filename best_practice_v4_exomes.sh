#!/bin/bash

# ******** Run parameters ********
BWA_DIR=/raid1/sequencing/apps/alignment/bwa-0.6.2
SAM_DIR=/raid1/sequencing/apps/post_processing/samtools-0.1.18
PICARD_DIR=/raid1/sequencing/apps/post_processing/picard-tools-1.78
GATK_DIR=/raid1/sequencing/apps/post_processing/GenomeAnalysisTK-2.1-11-g13c0244
REF_DIR=/raid1/sequencing/reference/gatk_bundle/gatk_bundle_04_Oct_2012

REF_FASTA=$REF_DIR/human_g1k_v37.fasta
# this only needs to be run once per reference .fasta
#$BWA_DIR/bwa index $REF_FASTA

DATA_DIR=/raid1/alex/sequencing/cll/data/exomes
SAMPLES=(`ls $DATA_DIR`)

NUM_CORES=12
RUN_JAVA="/usr/local/java/jdk1.6.0_22_x64/bin/java -Xmx16g"

CURRENT_DATE=`date +"%d_%b_%Y"`
TARGET_DIR=/raid1/alex/sequencing/cll/runs/exome_$CURRENT_DATE

# ******** Uncomment these lines to start fresh (for debugging) ********
rm -rf $TARGET_DIR
mkdir $TARGET_DIR
# ********************************

echo "********************"
echo "Aligning with BWA..."
echo "********************"
rm -rf $TARGET_DIR/alignment
mkdir $TARGET_DIR/alignment

for i in ${SAMPLES[*]}
do
	if [[ -d $DATA_DIR/$i ]]
	then
		echo "Aligning "$i
		
		cd $DATA_DIR/$i
		read1=(`ls *R1*.gz`)
		read2=(`ls *R2*.gz`)
		cd $TARGET_DIR
		mkdir alignment/$i
		mkdir alignment/$i/logs
		readGroup="@RG\tID:"$i"_exome_1\tSM:"$i"\tPL:ILLUMINA\tLB:exome"
		
		# build alignment indexes
		counter=0
		for j in ${read1[*]}
		do
			# now align the stuff
			# don't background it since it can use multiple cores inside
			$BWA_DIR/bwa aln \
				-t $NUM_CORES \
				$REF_FASTA \
				$DATA_DIR/$i/${read1[$counter]} \
				>$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair1.sai \
				2>$TARGET_DIR/alignment/$i/logs/bwa_aln_read1.log
			$BWA_DIR/bwa aln \
				-t $NUM_CORES \
				$REF_FASTA \
				$DATA_DIR/$i/${read2[$counter]} \
				>$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair2.sai \
				2>$TARGET_DIR/alignment/$i/logs/bwa_aln_read2.log
			counter=$(( $counter+1 ))
		done
		echo "Building, sorting "$i" .bam files"
		# build and sort the .bam files
		counter=0
		for j in ${read1[*]}
		do
			# Three things going on here at once: bwa sampe pipes .sam lines to
			# samtools view, converting them to .bam lines, and pipes those to
			# samtools sort, which finally creates the resulting .bam file
			
			# It looks like Jacob did some testing on running these in parallel:
			# "... this requires that we don't set the memory allocation for each job too high,
			# 2GB for the samtools sort step seems to work, other piped elements set their own....
			# the memory allocation here will blow up if there are more than about 16 total pairs
			# of fastq files for the genome"
			$BWA_DIR/bwa sampe \
				-r "$readGroup" \
				$REF_FASTA \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair1.sai \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair2.sai \
				$DATA_DIR/$i/${read1[$counter]} \
				$DATA_DIR/$i/${read2[$counter]} | \
			$SAM_DIR/samtools view \
				-buS \
				-t $REF_FASTA \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).bam
			$SAM_DIR/samtools sort \
				-m 20000000000 \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).bam \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).sorted.bam &
		done
		
		# wait for those last jobs to complete before continuing
		FAIL=0
		for job in `jobs -p`
		do
			wait $job || let "FAIL+=1"
		done
		if [[ $FAIL == 0 ]]
		then
			echo $i" completed."
		else
			echo "A job failed!"
			exit 1
		fi
		
		# TODO: When doing genome runs, I need to merge genome lanes into a single .bam
	fi
done

echo "******************"
echo "Marking Duplicates"
echo "******************"
rm -rf $TARGET_DIR/dedup
mkdir $TARGET_DIR/dedup
mkdir $TARGET_DIR/dedup/logs
for i in ${SAMPLES[*]}
do
	# TODO: the genomes run will have a different INPUT file here
	$RUN_JAVA $PICARD_DIR/MarkDuplicates.jar \
		INPUT=$TARGET_DIR/alignment/$i/lane1.sorted.bam \
		OUTPUT=$TARGET_DIR/dedup/$i.bam \
		METRICS=$TARGET_DIR/dedup/logs/$i.metrics.txt \
		>$TARGET_DIR/dedup/logs/$i.log \
		2>$TARGET_DIR/dedup/logs/$i.err.log &
done

# wait for jobs to complete before continuing
FAIL=0
for job in `jobs -p`
do
	wait $job || let "FAIL+=1"
done
if [[ $FAIL == 0 ]]
then
	echo $i" completed."
else
	echo "A job failed!"
	exit 1
fi