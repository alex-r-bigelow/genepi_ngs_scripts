#!/bin/bash
# Aligns, post-processes, and calls Paired-End Illumina Exome reads according to the
# GATK best practice version 4 (Best, does not use ReduceReads or VariantRecalibrator - applies hard filters).
# By default, UnifiedGenotyper is used for variant calling, but this can be changed by supplying a number as
# a parameter to this script (the number, if it exists, will be used as the --minPruning option for haplotypeCaller)

# You should never run this script directly; see example.sh for an example of how to launch this script; you should make
# a copy of example.sh and tweak it for every distinct run you do

# ******** Helper Functions ********

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
		echo "A job failed in the $PHASE_START step!"
		echo "A job failed in the $PHASE_START step!" >&2
		exit 1
	fi
}

# this allows the script to end early (when the PHASE_STOP parameter is set to something other than "")
finish ()
{
	echo "Finished jobs"
	if [ $PHASE_STOP != "" ]
	then
		echo "Stopped before the $PHASE_STOP step"
	fi
	exit 1
}

# These values should be set by a temporary script that launches this one:
# (I may add more crap if I add ANNOVAR to the mix)

# PHASE_START		"setup", "align", "build_bam", "sort_bam", "post_process", "call", "filter", "annotate"
# PHASE_STOP		"align", "build_bam", "sort_bam", "post_process", "call", "filter", "annotate", ""
#
# DATA_DIR			Directory. Should contain directories for each sample (named appropriately - these
#					names will be reused all the way past the call phase). Each subdirectory should
#					contain paired *.R1.*fastq.gz and *.R2.*fastq.gz (each lane should have exactly
#					the same name except for the R1 and R2)
# TARGET_DIR		The directory in which to store files, run everything. You can safely experiment
#					with MIN_PRUNING and/or FILTER_MODE on the same directory; otherwise you should only
#					restart on existing data if something failed
# MIN_PRUNING		blank or an integer. If blank, unifiedGenotyper will be used for variant calling,
#					otherwise haplotypeCaller will be used with --minPruning MIN_PRUNING
#					(a low value for MIN_PRUNING will be very slow)
# EXOME_OR_GENOME	"exome" or "genome"
#					will follow the appropriate parameters, depending on whether we're dealing with
#					exome/targeted or whole genome.
# EXOME_TARGETS		path to .bed file
#					If you're doing exome/targeted sequencing, you must provide a .bed file for
#					the targeted regions
# NUM_CORES			Number of processors on the machine running jobs
# MAX_MEM			max GB RAM to allow these jobs to consume
# FILTER_MODE		"ALL", "PURGED", or "PASS" - "ALL" will run annotations/summaries on all
#					variants, "PURGED" will run annotations on only on-target variants for
#					exome/targeted studies, and "PASS" will run annotations only on variants that
#					PASS all filters (except the DP Filter, which is applied by hand afterward)
# VST_SET_OP		String to give VST for combining individuals

# Paths to the main directory of each app (I allow multiple paths to VAAST programs as I've been using the
# bleeding-edge beta version and I sometimes have to experiment with different versions):
# JAVA
# BWA_DIR
# SAM_DIR
# PICARD_DIR
# GATK_DIR
# SNPEFF_DIR
# VAAST
# VAAST_vaast_converter
# VAAST_VST
# VAAST_vaast_sort_gff
# VAAST_VAT
# VAAST_QC			( for quality-check.pl, but there are some bad chars there )

# Paths to reference files:
# --- Reference genomes ----
# REF_FASTA
# VAAST_FASTA		VAAST usually doesn't cooperate well with reference genomes other than what they give you;
#					if you run buildVAASTreference.py, you might get it to cooperate with the same reference
#					genome that you use for the rest of the pipeline... but no guarantees here
# --- GATK bundle ---
# REF_DBSNP_129
# REF_DBSNP
# REF_MILLS
# REF_KGP
# REF_HAPMAP
# REF_OMNI
# ---- snpEff ----
# snpEff's files are a little tricky - for the main snpeff run, you need to download its
# main database via snpEff on the command line:
#
# -jar snpEff.jar download -v hg19
#
# It will store it in snpEff's directory structure. These other two you need to download
# manually (and you can put them wherever you like)
# GWAS_CAT
# DBNSFP
# ---- VAAST ----
# VAAST_BACKGROUND
# VAAST_FEATURES


# ******** Some initial variable verifying/twisting that always happens ********
# verify all our required parameters are valid
if	! ([ "$PHASE_START" == "setup" ] || \
		[ "$PHASE_START" == "align" ] || \
		[ "$PHASE_START" == "build_bam" ] || \
		[ "$PHASE_START" == "sort_bam" ] || \
		[ "$PHASE_START" == "dedup" ] || \
		[ "$PHASE_START" == "realign" ] || \
		[ "$PHASE_START" == "recalibrate" ] || \
		[ "$PHASE_START" == "call" ] || \
		[ "$PHASE_START" == "filter" ] || \
		[ "$PHASE_START" == "annotate" ]) || \
	! ([ "$PHASE_STOP" == "align" ] || \
		[ "$PHASE_STOP" == "build_bam" ] || \
		[ "$PHASE_STOP" == "sort_bam" ] || \
		[ "$PHASE_STOP" == "dedup" ] || \
		[ "$PHASE_STOP" == "realign" ] || \
		[ "$PHASE_STOP" == "recalibrate" ] || \
		[ "$PHASE_STOP" == "call" ] || \
		[ "$PHASE_STOP" == "filter" ] || \
		[ "$PHASE_STOP" == "annotate" ] || \
		[ "$PHASE_STOP" == "" ]) || \
	[ ! -d $DATA_DIR ] || \
	[ ! -d $TARGET_DIR ] || \
	! (([ "$EXOME_OR_GENOME" == "exome" ] && [ -e $EXOME_TARGETS ]) || [ "$EXOME_OR_GENOME" == "genome" ]) || \
	[ -z $NUM_CORES ] || \
	[ -z $MAX_MEM ] || \
	! ([ "$FILTER_MODE" == "ALL" ] || \
		[ "$FILTER_MODE" == "PURGED" ] || \
		[ "$FILTER_MODE" == "PASS" ]) || \
	[ -z $VST_SET_OP ] || \
	[ ! -e $JAVA/bin/java ] || \
	[ ! -e $BWA_DIR/bwa ] || \
	[ ! -e $SAM_DIR/samtools ] || \
	[ ! -e $PICARD_DIR/MarkDuplicates.jar ] || \
	[ ! -e $GATK_DIR/GenomeAnalysisTK.jar ] || \
	[ ! -e $SNPEFF_DIR/snpEff.jar ] || \
	[ ! -e $VAAST_VAAST ] || \
	[ ! -e $VAAST_vaast_converter ] || \
	[ ! -e $VAAST_VST ] || \
	[ ! -e $VAAST_vaast_sort_gff ] || \
	[ ! -e $VAAST_VAT ] || \
	[ ! -e $VAAST_QC ] || \
	[ ! -e $REF_FASTA ] || \
	[ ! -e $VAAST_FASTA ] || \
	[ ! -e $REF_DBSNP_129 ] || \
	[ ! -e $REF_DBSNP ] || \
	[ ! -e $REF_MILLS ] || \
	[ ! -e $REF_KGP ] || \
	[ ! -e $REF_HAPMAP ] || \
	[ ! -e $REF_OMNI ] || \
	[ ! -e $GWAS_CAT ] || \
	[ ! -e $DBNSFP ] || \
	[ ! -e $VAAST_BACKGROUND ] || \
	[ ! -e $VAAST_FEATURES ]
then
	echo "Something is wrong with your parameters!"
	echo "Something is wrong with your parameters!" >&2
	exit 1
fi

# get the list of our samples
TEMP=(`ls $DATA_DIR`)
SAMPLES=""
NUM_SAMPLES=0
for s in ${TEMP[*]}
do
	if [[ -d $DATA_DIR/$s ]]
	then
		SAMPLES="$SAMPLES $s"
		NUM_SAMPLES=$(( $NUM_SAMPLES+1 ))
	fi
done

# some simpler memory shortcuts
MAX_BYTE_MEM=$(( $MAX_MEM*1024*1024*1024 ))
MAX_SAMPLE_MEM=$(( $MAX_BYTE_MEM/$NUM_SAMPLES ))

# a java shortcut to use the appropriate amount of memory
RUN_JAVA="$JAVA/bin/java -Xmx$MAX_SAMPLE_MEM"

VCF_NAME=""
if [ -z $MIN_PRUNING ]
then
	VCF_NAME="unifiedGenotyper"
else
	VCF_NAME="haplotypeCaller$MIN_PRUNING"
fi

# silently test if we have enough memory... but fail early if we don't
$RUN_JAVA >/dev/null 2>/dev/null &
waitForJobs

#******** Actual jobs start here ********

if [ "$PHASE_START" == "setup" ]
then
	echo "setup..."
	
	# index the REF_FASTA for bwa, samtools if needed
	if	[ ! -e $REF_FASTA.amb ] || \
		[ ! -e $REF_FASTA.ann ] || \
		[ ! -e $REF_FASTA.bwt ] || \
		[ ! -e $REF_FASTA.pac ] || \
		[ ! -e $REF_FASTA.sa ]
	then
		$BWA_DIR/bwa index $REF_FASTA &
	fi
	if [ ! -e $REF_FASTA.fai ]
	then
		$SAM_DIR/samtools faidx $REF_FASTA &
	fi
	
	# TODO: download snpEff database for hg19 here instead of making the user do it...
	# this is actually really buggy! it tries to download in the wrong places, and you
	# have to tweak a line in snpEff.config to get it to work
	#if [ `ls snpEff_*_hg19.zip | wc -l` > 0 ]
	#then
	#	MYDIR=`pwd`
	#	cd $SNPEFF_DIR
	#	$RUN_JAVA -jar snpEff.jar download -v hg19 &
	#	cd $MYDIR
	#fi
	
	# start totally fresh
	rm -rf $TARGET_DIR &
	waitForJobs
	mkdir $TARGET_DIR
	
	export PHASE_START="align"
fi

if [ "$PHASE_STOP" == "align" ]
then
	finish
fi

if [ "$PHASE_START" == "align" ]
then
	echo "align..."
	rm -rf $TARGET_DIR/alignment
	mkdir $TARGET_DIR/alignment
	
	for i in ${SAMPLES[*]}
	do
		echo "...bwa aln "$i
		mkdir $TARGET_DIR/alignment/$i
		mkdir $TARGET_DIR/alignment/$i/logs
		
		# build alignment indexes
		counter=0
		read1=(`ls $DATA_DIR/$i/*R1*.gz`)
		read2=(`ls $DATA_DIR/$i/*R2*.gz`)
		for j in ${read1[*]}
		do
			# now align the stuff
			# don't run these jobs together since they can use multiple cores internally
			$BWA_DIR/bwa aln \
				-t $NUM_CORES \
				$REF_FASTA \
				${read1[$counter]} \
				>$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair1.sai \
				2>$TARGET_DIR/alignment/$i/logs/bwa_aln_read1.log &
			waitForJobs
			$BWA_DIR/bwa aln \
				-t $NUM_CORES \
				$REF_FASTA \
				${read2[$counter]} \
				>$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair2.sai \
				2>$TARGET_DIR/alignment/$i/logs/bwa_aln_read2.log &
			counter=$(( $counter+1 ))
			waitForJobs
		done
	done
	
	export PHASE_START="build_bam"
fi

if [ "$PHASE_STOP" == "build_bam" ]
then
	finish
fi

if [ "$PHASE_START" == "build_bam" ]
then
	echo "build_bam..."
	for i in ${SAMPLES[*]}
	do
		readGroup="@RG\tID:"$i"_$EXOME_OR_GENOME\_1\tSM:"$i"\tPL:ILLUMINA\tLB:$EXOME_OR_GENOME"
		counter=0
		read1=(`ls $DATA_DIR/$i/*R1*.gz`)
		read2=(`ls $DATA_DIR/$i/*R2*.gz`)
		for j in ${read1[*]}
		do
			echo "...bwa sampe | samtools view: $i, lane "$(( $counter+1 ))
			$BWA_DIR/bwa sampe \
				-r "$readGroup" \
				$REF_FASTA \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair1.sai \
				$TARGET_DIR/alignment/$i/lane$(( $counter+1 )).pair2.sai \
				${read1[$counter]} \
				${read2[$counter]} \
				2>$TARGET_DIR/alignment/$i/logs/lane$(( $counter+1 )).sampe.err.log \
				| \
			$SAM_DIR/samtools view \
				-bSh \
				-F 0x04 \
				-t $REF_FASTA.fai \
				-o $TARGET_DIR/alignment/$i/lane$(( $counter+1 )).unsorted.bam \
				- \
				>$TARGET_DIR/alignment/$i/logs/lane$(( $counter+1 )).view.log \
				2>$TARGET_DIR/alignment/$i/logs/lane$(( $counter+1 )).view.err.log &
			counter=$(( $counter+1 ))
		done
		waitForJobs
	done
	
	export PHASE_START="sort_bam"
fi

if [ "$PHASE_STOP" == "sort_bam" ]
then
	finish
fi

if [ "$PHASE_START" == "sort_bam" ]
then
	echo "sort_bam..."
	for i in ${SAMPLES[*]}
	do
		NUMLANES=`ls $DATA_DIR/$i/*R1*.gz | wc -l`
		BYTES_PER_LANE=$(( $MAX_BYTE_MEM / $NUMLANES ))
		for ((j=1; j<=$NUMLANES; j++))
		do
			echo "...samtools sort: $i, lane $j"
			$SAM_DIR/samtools sort \
				$TARGET_DIR/alignment/$i/lane$j.unsorted.bam \
				$TARGET_DIR/alignment/$i/lane$j \
				>$TARGET_DIR/alignment/$i/logs/lane$j.sort.log \
				2>$TARGET_DIR/alignment/$i/logs/lane$j.sort.err.log &
			waitForJobs
		done
	done
	
	for i in ${SAMPLES[*]}
	do
		NUMLANES=`ls $DATA_DIR/$i/*R1*.gz | wc -l`
		if [ $NUMLANES == 1 ]
		then
			mv $TARGET_DIR/alignment/$i/lane1.bam $TARGET_DIR/alignment/$i/merged.bam
			continue
		fi
		echo "...samtools merge: $i"
		ALLLANES=""
		for ((j=1; j<=$NUMLANES; j++))
		do
			ALLLANES="$ALLLANES $TARGET_DIR/alignment/$i/lane$j.bam"
		done
		$SAM_DIR/samtools merge \
			$TARGET_DIR/alignment/$i/merged.bam \
			$ALLLANES \
			>$TARGET_DIR/alignment/$i/logs/merge.log \
			2>$TARGET_DIR/alignment/$i/logs/merge.err.log &
		waitForJobs
	done
	
	export PHASE_START="dedup"
fi

if [ "$PHASE_STOP" == "dedup" ]
then
	finish
fi

if [ "$PHASE_START" == "dedup" ]
then
	echo "dedup..."
	rm -rf $TARGET_DIR/dedup
	mkdir $TARGET_DIR/dedup
	mkdir $TARGET_DIR/dedup/logs
	for i in ${SAMPLES[*]}
	do
		echo "...Picard MarkDuplicates: "$i
		$RUN_JAVA -jar $PICARD_DIR/MarkDuplicates.jar \
			INPUT=$TARGET_DIR/alignment/$i/merged.bam \
			OUTPUT=$TARGET_DIR/dedup/$i.bam \
			METRICS_FILE=$TARGET_DIR/dedup/logs/$i.metrics.txt \
			AS=true \
			VALIDATION_STRINGENCY=LENIENT \
			>$TARGET_DIR/dedup/logs/$i.log \
			2>$TARGET_DIR/dedup/logs/$i.err.log &
	done
	waitForJobs
	
	for i in ${SAMPLES[*]}
	do
		echo "...samtools index: "$i
		$SAM_DIR/samtools index \
			$TARGET_DIR/dedup/$i.bam \
			$TARGET_DIR/dedup/$i.bai &
	done
	waitForJobs
	
	export PHASE_START="realign"
fi

if [ "$PHASE_STOP" == "realign" ]
then
	finish
fi

if [ "$PHASE_START" == "realign" ]
then
	echo "realign..."
	rm -rf $TARGET_DIR/realignment
	mkdir $TARGET_DIR/realignment
	mkdir $TARGET_DIR/realignment/logs
	mkdir $TARGET_DIR/realignment/stats
	for i in ${SAMPLES[*]}
	do
		echo "...GATK RealignerTargetCreator: "$i
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T RealignerTargetCreator \
			-I $TARGET_DIR/dedup/$i.bam \
			-o $TARGET_DIR/realignment/stats/$i.intervals \
			-R $REF_FASTA \
			--known $REF_MILLS \
			--known $REF_KGP \
			>$TARGET_DIR/realignment/logs/$i.realignerTargetCreator.log \
			2>$TARGET_DIR/realignment/logs/$i.realignerTargetCreator.err.log &
	done
	waitForJobs
	
	for i in ${SAMPLES[*]}
	do
		echo "...GATK IndelRealigner: "$i
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-I $TARGET_DIR/dedup/$i.bam \
			-o $TARGET_DIR/realignment/$i.bam \
			-R $REF_FASTA \
			-targetIntervals $TARGET_DIR/realignment/stats/$i.intervals \
			-known $REF_MILLS \
			-known $REF_KGP \
			>$TARGET_DIR/realignment/logs/$i.indelRealigner.log \
			2>$TARGET_DIR/realignment/logs/$i.indelRealigner.err.log &
	done
	waitForJobs

	for i in ${SAMPLES[*]}
	do
		echo "...samtools index: "$i
		$SAM_DIR/samtools index \
			$TARGET_DIR/realignment/$i.bam \
			$TARGET_DIR/realignment/$i.bai &
	done
	waitForJobs
	
	export PHASE_START="recalibrate"
fi

if [ "$PHASE_STOP" == "recalibrate" ]
then
	finish
fi

if [ "$PHASE_START" == "recalibrate" ]
then
	rm -rf $TARGET_DIR/recalibration
	mkdir $TARGET_DIR/recalibration
	mkdir $TARGET_DIR/recalibration/logs
	mkdir $TARGET_DIR/recalibration/stats
	
	echo "recalibrate..."
	for i in ${SAMPLES[*]}
	do
		echo "...GATK BaseRecalibrator: "$i
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-I $TARGET_DIR/realignment/$i.bam \
			-o $TARGET_DIR/recalibration/stats/$i.grp \
			-R $REF_FASTA \
			-knownSites $REF_MILLS \
			-knownSites $REF_KGP \
			-knownSites $REF_DBSNP \
			>$TARGET_DIR/recalibration/logs/$i.baseRecalibrator.log \
			2>$TARGET_DIR/recalibration/logs/$i.baseRecalibrator.err.log &
	done
	waitForJobs
	
	for i in ${SAMPLES[*]}
	do
		echo "...GATK PrintReads (update scores): "$i
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T PrintReads \
			-I $TARGET_DIR/realignment/$i.bam \
			-o $TARGET_DIR/recalibration/$i.bam \
			-R $REF_FASTA \
			-BQSR $TARGET_DIR/recalibration/stats/$i.grp \
			>$TARGET_DIR/recalibration/logs/$i.updateScores.log \
			2>$TARGET_DIR/recalibration/logs/$i.updateScores.err.log &
	done
	waitForJobs
	
	for i in ${SAMPLES[*]}
	do
		echo "...samtools index: "$i
		$SAM_DIR/samtools index \
			$TARGET_DIR/recalibration/$i.bam \
			$TARGET_DIR/recalibration/$i.bai &
	done
	waitForJobs
	
	rm -rf $TARGET_DIR/calls
	mkdir $TARGET_DIR/calls
	mkdir $TARGET_DIR/calls/logs
	
	rm -rf $TARGET_DIR/annotation
	mkdir $TARGET_DIR/annotation
	mkdir $TARGET_DIR/annotation/vaast
	mkdir $TARGET_DIR/annotation/snpeff
	mkdir $TARGET_DIR/annotation/snpeff/logs
	
	TEMP=""
	for i in ${SAMPLES[*]}
	do
		TEMP="$TEMP -I $TARGET_DIR/recalibration/$i.bam"
	done
	echo "...GATK PrintReads (merge):"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T PrintReads \
		$TEMP \
		-o $TARGET_DIR/calls/all.bam \
		-R $REF_FASTA \
		>$TARGET_DIR/calls/logs/all.merge.log \
		2>$TARGET_DIR/calls/logs/all.merge.err.log &
	waitForJobs
	
	export PHASE_START="call"
fi

if [ "$PHASE_STOP" == "call" ]
then
	finish
fi

if [ "$PHASE_START" == "call" ]
then
	echo "call..."
	
	if [ -z $MIN_PRUNING ]
	then
		echo "...GATK UnifiedGenotyper (SNPs)"
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T UnifiedGenotyper \
			-I $TARGET_DIR/calls/all.bam \
			-o $TARGET_DIR/calls/snps.$VCF_NAME.raw.vcf \
			-R $REF_FASTA \
			--dbsnp $REF_DBSNP \
			--genotype_likelihoods_model SNP \
			>$TARGET_DIR/calls/logs/snps.$VCF_NAME.raw.log \
			2>$TARGET_DIR/calls/logs/snps.$VCF_NAME.raw.err.log &
		echo "...GATK UnifiedGenotyper (INDELs)"
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T UnifiedGenotyper \
			-I $TARGET_DIR/calls/all.bam \
			-o $TARGET_DIR/calls/indels.$VCF_NAME.raw.vcf \
			-R $REF_FASTA \
			--dbsnp $REF_DBSNP \
			--genotype_likelihoods_model INDEL \
			>$TARGET_DIR/calls/logs/indels.$VCF_NAME.raw.log \
			2>$TARGET_DIR/calls/logs/indels.$VCF_NAME.raw.err.log &
		waitForJobs
	else
		echo "...GATK HaplotypeCaller --minPruning $MIN_PRUNING"
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-I $TARGET_DIR/calls/all.bam \
			-o $TARGET_DIR/calls/$VCF_NAME.raw.vcf \
			-R $REF_FASTA \
			--dbsnp $REF_DBSNP \
			--minPruning $MIN_PRUNING \
			>$TARGET_DIR/calls/logs/$VCF_NAME.log \
			2>$TARGET_DIR/calls/logs/$VCF_NAME.err.log &
		waitForJobs
		
		echo "...GATK SelectVariants (separate SNPs, INDELs)"
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T SelectVariants \
			--variant $TARGET_DIR/calls/$VCF_NAME.raw.vcf \
			-o $TARGET_DIR/calls/indels.$VCF_NAME.raw.vcf \
			-R $REF_FASTA \
			-selectType INDEL \
			>$TARGET_DIR/calls/logs/indels.$VCF_NAME.selectVariants.log \
			2>$TARGET_DIR/calls/logs/indels.$VCF_NAME.selectVariants.err.log &
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T SelectVariants \
			--variant $TARGET_DIR/calls/$VCF_NAME.raw.vcf \
			-o $TARGET_DIR/calls/snps.$VCF_NAME.raw.vcf \
			-R $REF_FASTA \
			-selectType SNP \
			>$TARGET_DIR/calls/logs/snps.$VCF_NAME.selectVariants.log \
			2>$TARGET_DIR/calls/logs/snps.$VCF_NAME.selectVariants.err.log &
		waitForJobs
	fi
	
	export PHASE_START="filter"
fi

if [ "$PHASE_STOP" == "filter" ]
then
	finish
fi

if [ "$PHASE_START" == "filter" ]
then
	echo "filter..."
	
	INBREEDING_PARAMETER=""
	if (( ${#SAMPLES[@]} >= 10 ))
	then
		INBREEDING_PARAMETER="--filterExpression \"InbreedingCoeff < -0.8\" --filterName \"InbreedingCoeff Filter\""
	fi
	echo "...GATK VariantFiltration (SNPs)"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		--variant $TARGET_DIR/calls/snps.$VCF_NAME.raw.vcf \
		-o $TARGET_DIR/calls/snps.$VCF_NAME.filtered.vcf \
		-R $REF_FASTA \
		--filterExpression "QD < 2.0" \
		--filterName "QD Filter" \
		--filterExpression "MQ < 40.0" \
		--filterName "MQ Filter" \
		--filterExpression "FS > 60.0" \
		--filterName "FS Filter" \
		--filterExpression "HaplotypeScore > 13.0" \
		--filterName "HaplotypeScore Filter" \
		--filterExpression "MQRankSum < -12.5" \
		--filterName "MQRankSum Filter" \
		--filterExpression "ReadPosRankSum < -8.0" \
		--filterName "ReadPosRankSum Filter" \
		>$TARGET_DIR/calls/logs/snps.$VCF_NAME.variantFiltration.log \
		2>$TARGET_DIR/calls/logs/snps.$VCF_NAME.variantFiltration.err.log &
	echo "...GATK VariantFiltration (INDELs)"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		--variant $TARGET_DIR/calls/indels.$VCF_NAME.raw.vcf \
		-o $TARGET_DIR/calls/indels.$VCF_NAME.filtered.vcf \
		-R $REF_FASTA \
		--filterExpression "QD < 2.0" \
		--filterName "QD Filter" \
		--filterExpression "ReadPosRankSum < -20.0" \
		--filterName "ReadPosRankSum Filter" \
		--filterExpression "FS > 200.0" \
		--filterName "FS Filter" \
		$INBREEDING_PARAMETER \
		>$TARGET_DIR/calls/logs/indels.$VCF_NAME.variantFiltration.log \
		2>$TARGET_DIR/calls/logs/indels.$VCF_NAME.variantFiltration.err.log &
	waitForJobs
	
	echo "...GATK CombineVariants (recombine SNPs and INDELs)"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T CombineVariants \
		--variant $TARGET_DIR/calls/snps.$VCF_NAME.filtered.vcf \
		--variant $TARGET_DIR/calls/indels.$VCF_NAME.filtered.vcf \
		-o $TARGET_DIR/calls/all.$VCF_NAME.filtered.vcf \
		-R $REF_FASTA \
		>$TARGET_DIR/calls/logs/all.$VCF_NAME.combineVariants.log \
		2>$TARGET_DIR/calls/logs/all.$VCF_NAME.combineVariants.err.log &
	waitForJobs
	
	if [ "$EXOME_OR_GENOME" == "exome" ]
	then
		echo "...GATK SelectVariants (remove off-target variants)"
		$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
			-T SelectVariants \
			--variant $TARGET_DIR/calls/all.$VCF_NAME.filtered.vcf \
			-o $TARGET_DIR/calls/purged.$VCF_NAME.filtered.vcf \
			-R $REF_FASTA \
			-L $EXOME_TARGETS \
			>$TARGET_DIR/calls/logs/purged.$VCF_NAME.selectVariants.log \
			2>$TARGET_DIR/calls/logs/purged.$VCF_NAME.selectVariants.err.log &
		waitForJobs
	else
		cp $TARGET_DIR/calls/all.$VCF_NAME.filtered.vcf $TARGET_DIR/calls/purged.$VCF_NAME.filtered.vcf
	fi
	
	echo "...GATK SelectVariants (physically remove variants that didn't PASS)"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T SelectVariants \
		--variant $TARGET_DIR/calls/purged.$VCF_NAME.filtered.vcf \
		-o $TARGET_DIR/calls/pass.$VCF_NAME.filtered.vcf \
		-R $REF_FASTA \
		--excludeFiltered \
		>$TARGET_DIR/calls/logs/pass.$VCF_NAME.selectVariants.log \
		2>$TARGET_DIR/calls/logs/pass.$VCF_NAME.selectVariants.err.log &
	waitForJobs
	
	echo "...dpFilter.py (manual, std-dev-based filter that remains in the .vcf)"
	python ${0%best_practice_v4.sh}dpFilter.py \
		--stddev 5 \
		--in $TARGET_DIR/calls/all.$VCF_NAME.filtered.vcf \
		--out $TARGET_DIR/calls/ALL.$VCF_NAME.vcf \
		>$TARGET_DIR/calls/logs/all.$VCF_NAME.dpFilter.log \
		2>$TARGET_DIR/calls/logs/all.$VCF_NAME.dpFilter.err.log &
	python ${0%best_practice_v4.sh}dpFilter.py \
		--stddev 5 \
		--in $TARGET_DIR/calls/purged.$VCF_NAME.filtered.vcf \
		--out $TARGET_DIR/calls/PURGED.$VCF_NAME.vcf \
		>$TARGET_DIR/calls/logs/purged.$VCF_NAME.dpFilter.log \
		2>$TARGET_DIR/calls/logs/purged.$VCF_NAME.dpFilter.err.log &
	python ${0%best_practice_v4.sh}dpFilter.py \
		--stddev 5 \
		--in $TARGET_DIR/calls/pass.$VCF_NAME.filtered.vcf \
		--out $TARGET_DIR/calls/PASS.$VCF_NAME.vcf \
		>$TARGET_DIR/calls/logs/pass.$VCF_NAME.dpFilter.log \
		2>$TARGET_DIR/calls/logs/pass.$VCF_NAME.dpFilter.err.log &
	waitForJobs
	export PHASE_START="annotate"
fi

# At this point, the .vcf files with CAPS in their names are ready for analysis elsewhere... FILTER_MODE will decide which one gets annotated
VCF_NAME=$FILTER_MODE.$VCF_NAME

if [ "$PHASE_STOP" == "annotate" ]
then
	finish
fi

if [ "$PHASE_START" == "annotate" ]
then
	echo "annotate..."
	
	
	echo "...VAAST"
	rm -rf $TARGET_DIR/annotation/vaast/$VCF_NAME
	mkdir $TARGET_DIR/annotation/vaast/$VCF_NAME
	mkdir $TARGET_DIR/annotation/vaast/$VCF_NAME/logs
	
	echo "......removeAlternateContigs"
	python ${0%best_practice_v4.sh}removeAlternateContigs.py \
		--in $TARGET_DIR/calls/$VCF_NAME.vcf \
		--out $TARGET_DIR/annotation/vaast/$VCF_NAME/basicContigs.vcf \
		>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/removeAlternateContigs.log \
		2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/removeAlternateContigs.err.log &
	waitForJobs
	
	echo "......vaast_converter"
	mkdir $TARGET_DIR/annotation/vaast/$VCF_NAME/pre
	$VAAST_vaast_converter \
		--build hg19 \
		--format VCF \
		--path $TARGET_DIR/annotation/vaast/$VCF_NAME/pre/ \
		$TARGET_DIR/annotation/vaast/$VCF_NAME/basicContigs.vcf \
		>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/vaast_converter.log \
		2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/vaast_converter.err.log &
	waitForJobs
	
	echo "......vaast_sort_gff3"
	for i in $TARGET_DIR/annotation/vaast/$VCF_NAME/pre/*.gvf
	do
		j=${i##*/}
		$VAAST_vaast_sort_gff \
			--in_place \
			--perl_sort \
			--no_backup \
			$i \
			>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/${j%.gvf}.pre.vaast_sort_gff.log \
			2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/${j%.gvf}.pre.vaast_sort_gff.err.log &
	done
	waitForJobs
	
	echo "......VAT"
	rm -rf $TARGET_DIR/annotation/vaast/$VCF_NAME/post
	mkdir $TARGET_DIR/annotation/vaast/$VCF_NAME/post
	for i in $TARGET_DIR/annotation/vaast/$VCF_NAME/pre/*.gvf
	do
		j=${i##*/}
		$VAAST_VAT \
			--build hg19 \
			--fasta $VAAST_FASTA \
			--features $VAAST_FEATURES \
			$i \
			>$TARGET_DIR/annotation/vaast/$VCF_NAME/post/$j \
			2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/${j%.gvf}.VAT.log &
	done
	waitForJobs
	
	echo "......vaast_sort_gff3"
	for i in $TARGET_DIR/annotation/vaast/$VCF_NAME/post/*.gvf
	do
		j=${i##*/}
		$VAAST_vaast_sort_gff \
			--in_place \
			--perl_sort \
			--no_backup \
			$i \
			>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/${j%.gvf}.post.vaast_sort_gff.log \
			2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/${j%.gvf}.post.vaast_sort_gff.err.log &
	done
	waitForJobs
	
	echo "......VST"
	echo $VST_SET_OP
	$VAAST_VST \
		--ops $VST_SET_OP \
		$TARGET_DIR/annotation/vaast/$VCF_NAME/post/*.gvf \
		>$TARGET_DIR/annotation/vaast/$VCF_NAME/$VCF_NAME.cdr \
		2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/VST.log &
	waitForJobs
	
	echo "......quality-check.pl"
	$VAAST_QC \
		-sim 100000 \
		$VAAST_BACKGROUND \
		$TARGET_DIR/annotation/vaast/$VCF_NAME/$VCF_NAME.cdr \
		>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/quality-check.log \
		2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/quality-check.err.log &
	waitForJobs
	
	echo "......VAAST"
	$VAAST \
		--mode lrt \
		--outfile $TARGET_DIR/annotation/vaast/$VCF_NAME/$VCF_NAME \
		--indel \
		$VAAST_FEATURES \
		$VAAST_BACKGROUND \
		$TARGET_DIR/annotation/vaast/$VCF_NAME/$VCF_NAME.cdr \
		>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/VAAST.log \
		2>$TARGET_DIR/annotation/vaast/$VCF_NAME/logs/VAAST.err.log &
	waitForJobs
	
	echo "...snpeff"
	
	echo "......normal"
	$RUN_JAVA -jar $SNPEFF_DIR/snpEff.jar \
		-c $SNPEFF_DIR/snpEff.config \
		-s $TARGET_DIR/annotation/snpeff/$VCF_NAME.summary.html \
		hg19 \
		$TARGET_DIR/calls/$VCF_NAME.vcf \
		>$TARGET_DIR/annotation/snpeff/$VCF_NAME.predbnsfp.vcf \
		2>$TARGET_DIR/annotation/snpeff/logs/$VCF_NAME.predbnsfp.err.log &
	waitForJobs
	
	echo "......dbnsfp"
	$RUN_JAVA -jar $SNPEFF_DIR/SnpSift.jar \
		dbnsfp \
		-a \
		$DBNSFP \
		$TARGET_DIR/annotation/snpeff/$VCF_NAME.predbnsfp.vcf \
		>$TARGET_DIR/annotation/snpeff/$VCF_NAME.pregwascat.vcf \
		2>$TARGET_DIR/annotation/snpeff/logs/$VCF_NAME.pregwascat.err.log &
	waitForJobs
	
	echo "......gwasCat"
	$RUN_JAVA -jar $SNPEFF_DIR/SnpSift.jar \
		gwasCat \
		$GWAS_CAT \
		$TARGET_DIR/annotation/snpeff/$VCF_NAME.pregwascat.vcf \
		>$TARGET_DIR/annotation/snpeff/$VCF_NAME.snpeff.vcf \
		2>$TARGET_DIR/annotation/snpeff/logs/$VCF_NAME.snpeff.err.log &
	waitForJobs
fi

echo "Done"