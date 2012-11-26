#!/bin/bash
# Aligns, post-processes, and calls Paired-End Illumina Whole Genome reads according to the
# GATK best practice version 4 (Best, does not use ReduceReads or VariantRecalibrator - applies hard filters).
# By default, UnifiedGenotyper is used for variant calling, but this can be changed by supplying a number as
# a parameter to this script (the number, if it exists, will be used as the --minPruning option for haplotypeCaller)

# Redirect DATA_DIR and TARGET_DIR for another project (currently works for the cll project)


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
BWA_DIR=/raid1/sequencing/apps/alignment/bwa-0.6.2
SAM_DIR=/raid1/sequencing/apps/post_processing/samtools-0.1.18
PICARD_DIR=/raid1/sequencing/apps/post_processing/picard-tools-1.78
GATK_DIR=/raid1/sequencing/apps/post_processing/GenomeAnalysisTK-2.1-11-g13c0244
REF_DIR=/raid1/sequencing/reference/gatk_bundle/gatk_bundle_04_Oct_2012

REF_DBSNP_129=$REF_DIR/dbsnp_135.b37.excluding_sites_after_129.vcf
REF_DBSNP=$REF_DIR/dbsnp_135.b37.vcf
REF_MILLS=$REF_DIR/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
REF_KGP=$REF_DIR/1000G_phase1.indels.b37.vcf
REF_HAPMAP=$REF_DIR/hapmap_3.3.b37.vcf
REF_OMNI=$REF_DIR/1000G_omni2.5.b37.sites.vcf
REF_FASTA=$REF_DIR/human_g1k_v37.fasta
# <<<< these only need to be run once per reference .fasta >>>>
#$BWA_DIR/bwa index $REF_FASTA &
#$SAM_DIR/samtools faidx $REF_FASTA &
#waitForJobs

DATA_DIR=/raid1/alex/sequencing/cll/data/genomes
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

# <<<< adjust these two values to match the number of cores, GB of ram on the machine >>>>
NUM_CORES=12
MAX_MEM=64

MAX_BYTE_MEM=$(( $MAX_MEM*1024*1024*1024 ))
MAX_SAMPLE_MEM=$(( $MAX_BYTE_MEM/$NUM_SAMPLES ))
RUN_JAVA="/usr/local/java/jdk1.6.0_22_x64/bin/java -Xmx$MAX_SAMPLE_MEM"

$RUN_JAVA >/dev/null 2>/dev/null &
waitForJobs
# silently test if we have enough memory... but fail early if we don't

CURRENT_DATE=`date +"%d_%b_%Y"`
# <<<< uncomment this to add to an existing day's run >>>>
#CURRENT_DATE=22_Oct_2012
TARGET_DIR=/raid1/alex/sequencing/cll/runs/genome_$CURRENT_DATE

VCF_NAME=""
if [ -z $1 ]
then
	VCF_NAME="unifiedGenotyper"
else
	VCF_NAME="haplotypeCaller$1"
fi

# I use this if statement to comment stuff out
if [[ 0 == 1 ]]
then
# ******** Actual jobs start here ********
rm -rf $TARGET_DIR
mkdir $TARGET_DIR

echo "********************"
echo "Aligning with BWA..."
echo "********************"
rm -rf $TARGET_DIR/alignment
mkdir $TARGET_DIR/alignment
NUMLANES=0
for i in ${SAMPLES[*]}
do
	echo "Aligning "$i
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
		NUMLANES=$(( $NUMLANES+1 ))
		waitForJobs
	done
done

fi

echo "Building, per-lane .bam files"
for i in ${SAMPLES[*]}
do
	readGroup="@RG\tID:"$i"_genome_1\tSM:"$i"\tPL:ILLUMINA\tLB:genome"
	
	counter=0
	read1=(`ls $DATA_DIR/$i/*R1*.gz`)
	read2=(`ls $DATA_DIR/$i/*R2*.gz`)
	for j in ${read1[*]}
	do
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
	done
	waitForJobs
done

BYTES_PER_LANE=$(( $MAX_BYTE_MEM/$NUMLANES ))
echo "Sorting per-lane .bam files"
for i in ${SAMPLES[*]}
do
	NUMLANES=`ls $DATA_DIR/$i/*R1*.gz | wc -l`
	for ((j=1; j<=$NUMLANES; j++))
	do
		$SAM_DIR/samtools sort \
			-m 2000000000 \
			$TARGET_DIR/alignment/$i/lane$j.unsorted.bam \
			$TARGET_DIR/alignment/$i/lane$j \
			>$TARGET_DIR/alignment/$i/logs/lane$j.sort.log \
			2>$TARGET_DIR/alignment/$i/logs/lane$j.sort.err.log &
	done
done
waitForJobs

echo "Bundling per-sample .bam files"
for i in ${SAMPLES[*]}
do
	NUMLANES=`ls $DATA_DIR/$i/*R1*.gz | wc -l`
	if [ $NUMLANES == 1 ]
	then
		mv $TARGET_DIR/alignment/$i/lane1.bam $TARGET_DIR/alignment/$i/merged.bam
		continue
	fi
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
done
waitForJobs

echo "******************"
echo "Marking Duplicates"
echo "******************"
rm -rf $TARGET_DIR/dedup
mkdir $TARGET_DIR/dedup
mkdir $TARGET_DIR/dedup/logs
for i in ${SAMPLES[*]}
do
	# TODO: the genomes run will have a different INPUT file here
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
echo "Indexing..."
for i in ${SAMPLES[*]}
do
	$SAM_DIR/samtools index \
		$TARGET_DIR/dedup/$i.bam \
		$TARGET_DIR/dedup/$i.bai &
done
waitForJobs

echo "*************************************"
echo "Local realignment around small INDELS"
echo "*************************************"
rm -rf $TARGET_DIR/realignment
mkdir $TARGET_DIR/realignment
mkdir $TARGET_DIR/realignment/logs
mkdir $TARGET_DIR/realignment/stats
echo "Generating targets"
for i in ${SAMPLES[*]}
do
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

echo "Realigning"
for i in ${SAMPLES[*]}
do
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

echo "Indexing..."
for i in ${SAMPLES[*]}
do
	$SAM_DIR/samtools index \
		$TARGET_DIR/realignment/$i.bam \
		$TARGET_DIR/realignment/$i.bai &
done
waitForJobs

echo "**************************"
echo "Base Quality Recalibration"
echo "**************************"
rm -rf $TARGET_DIR/recalibration
mkdir $TARGET_DIR/recalibration
mkdir $TARGET_DIR/recalibration/logs
mkdir $TARGET_DIR/recalibration/stats

echo "Calculating Stats"
for i in ${SAMPLES[*]}
do
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-I $TARGET_DIR/realignment/$i.bam \
		-o $TARGET_DIR/recalibration/stats/$i.grp \
		-R $REF_FASTA \
		--knownSites $REF_MILLS \
		--knownSites $REF_KGP \
		--knownSites $REF_DBSNP \
		>$TARGET_DIR/recalibration/logs/$i.baseRecalibrator.log \
		2>$TARGET_DIR/recalibration/logs/$i.baseRecalibrator.err.log &
done
waitForJobs

echo "Adjusting scores"
for i in ${SAMPLES[*]}
do
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

echo "Indexing..."
for i in ${SAMPLES[*]}
do
	$SAM_DIR/samtools index \
		$TARGET_DIR/recalibration/$i.bam \
		$TARGET_DIR/recalibration/$i.bai &
done
waitForJobs

echo "***************"
echo "Variant Calling"
echo "***************"
rm -rf $TARGET_DIR/calls
mkdir $TARGET_DIR/calls
mkdir $TARGET_DIR/calls/logs

echo "Merging samples"
TEMP=""
for i in ${SAMPLES[*]}
do
	TEMP="$TEMP -I $TARGET_DIR/recalibration/$i.bam"
done
$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
	-T PrintReads \
	$TEMP \
	-o $TARGET_DIR/calls/all.bam \
	-R $REF_FASTA \
	>$TARGET_DIR/calls/logs/all.merge.log \
	2>$TARGET_DIR/calls/logs/all.merge.err.log &
waitForJobs

echo "Calling Variants using:"
if [ -z $1 ]
then
	echo "UnifiedGenotyper"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T UnifiedGenotyper \
		-I $TARGET_DIR/calls/all.bam \
		-o $TARGET_DIR/calls/snps.$VCF_NAME.raw.vcf \
		-R $REF_FASTA \
		--dbsnp $REF_DBSNP \
		--genotype_likelihoods_model SNP \
		>$TARGET_DIR/calls/logs/snps.$VCF_NAME.raw.log \
		2>$TARGET_DIR/calls/logs/snps.$VCF_NAME.raw.err.log &
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
	echo "HaplotypeCaller --minPruning $1"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-I $TARGET_DIR/calls/all.bam \
		-o $TARGET_DIR/calls/all.$VCF_NAME.raw.vcf \
		-R $REF_FASTA \
		--dbsnp $REF_DBSNP \
		--minPruning $1 \
		>$TARGET_DIR/calls/logs/all.$VCF_NAME.log \
		2>$TARGET_DIR/calls/logs/all.$VCF_NAME.err.log &
	waitForJobs
	
	echo "Separating SNPs, INDELs for filtering"
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T SelectVariants \
		--variant $TARGET_DIR/calls/all.$VCF_NAME.raw.vcf \
		-o $TARGET_DIR/calls/indels.$VCF_NAME.raw.vcf \
		-R $REF_FASTA \
		-selectType INDEL \
		>$TARGET_DIR/calls/logs/indels.$VCF_NAME.selectVariants.log \
		2>$TARGET_DIR/calls/logs/indels.$VCF_NAME.selectVariants.err.log &
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T SelectVariants \
		--variant $TARGET_DIR/calls/all.$VCF_NAME.raw.vcf \
		-o $TARGET_DIR/calls/snps.$VCF_NAME.raw.vcf \
		-R $REF_FASTA \
		-selectType SNP \
		>$TARGET_DIR/calls/logs/snps.$VCF_NAME.selectVariants.log \
		2>$TARGET_DIR/calls/logs/snps.$VCF_NAME.selectVariants.err.log &
	waitForJobs
fi

echo "Applying hard filters"
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
	>$TARGET_DIR/calls/logs/indels.$VCF_NAME.variantFiltration.log \
	2>$TARGET_DIR/calls/logs/indels.$VCF_NAME.variantFiltration.err.log &
waitForJobs

echo "Recombining filtered SNPs and INDELs"
$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
	-T CombineVariants \
	--variant $TARGET_DIR/calls/snps.$VCF_NAME.filtered.vcf \
	--variant $TARGET_DIR/calls/indels.$VCF_NAME.filtered.vcf \
	-o $TARGET_DIR/calls/all.$VCF_NAME.filtered.vcf \
	-R $REF_FASTA \
	>$TARGET_DIR/calls/logs/all.$VCF_NAME.combineVariants.log \
	2>$TARGET_DIR/calls/logs/all.$VCF_NAME.combineVariants.err.log &
waitForJobs

echo "Applying DP Filter by hand"
python dpFilter.py \
	--stddev 5 \
	--in $TARGET_DIR/calls/all.$VCF_NAME.filtered.vcf \
	--out $TARGET_DIR/calls/all.$VCF_NAME.finished.vcf \
	>$TARGET_DIR/calls/logs/all.$VCF_NAME.dpFilter.log \
	2>$TARGET_DIR/calls/logs/all.$VCF_NAME.dpFilter.err.log &
waitForJobs