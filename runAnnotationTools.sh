#!/bin/bash
# To be run after one of the best_practice_v4 pipelines; will add an additional directory containing various
# annotations. The script will run VAAST, calcStats.py, SnpEff, and ANNOVAR

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
# <<<< edit these values (and best_paths.txt) to point to the appropriate programs, reference files >>>>
while read p; do
	if [[ $p == *vaast_converter* ]]
	then
		VAAST_CONVERTER=$p
		echo "Using: "$VAAST_CONVERTER
	fi
	if [[ $p == *VST ]]
	then
		VST=$p
		echo "Using: "$VST
	fi
	if [[ $p == *VAAST ]]
	then
		VAAST=$p
		echo "Using: "$VAAST
	fi
	if [[ $p == *vaast_sort_gff ]]
	then
		VAAST_SORT_GFF=$p
		echo "Using: "$VAAST_SORT_GFF
	fi
	if [[ $p == *VAT ]]
	then
		VAT=$p
		echo "Using: "$VAT
	fi
	if [[ $p == *quality-check.pl ]]
	then
		QUALITY_CHECK_PL=$p
		echo "Using: "$QUALITY_CHECK_PL
	fi
done < best_paths.txt
if [ -z $QUALITY_CHECK_PL ] || [ -z $VAT ] || [ -z $VAAST_SORT_GFF ] || [ -z $VAAST ] || [ -z $VST ] || [ -z $VAAST_CONVERTER ]
then
	echo "Missing VAAST path(s)"
	exit 1
fi
GATK_DIR=/raid1/sequencing/apps/post_processing/GenomeAnalysisTK-2.1-11-g13c0244
SNPEFF_DIR=/raid1/sequencing/apps/annotation/snpEff_3_1
ANNOVAR_DIR=/raid1/sequencing/apps/annotation/annovar

KGP_DIR=/raid1/sequencing/reference/background/KGP/
KGP_DATA_DIR=$KGP_DIR/compressed_vcfs
KGP_POP_DIR=$KGP_DIR/populationLists

REF_GENE=/raid1/sequencing/reference/vaast/refGene_hg19.gff3

GWAS_CAT=/raid1/sequencing/reference/other/gwascatalog.txt
DBNSFP=/raid1/sequencing/reference/other/dbNSFP2.0b3.txt

VAAST_BACKGROUND=/raid1/sequencing/reference/background/vaast/1KG_refGene_Dec2011_CGDiv_NHLBI_NoCall.cdr

REF_DIR=/raid1/sequencing/reference/gatk_bundle/gatk_bundle_04_Oct_2012
REF_FASTA=$REF_DIR/human_g1k_v37.fasta
REF_FASTA2=$REF_DIR/human_g1k_v37.vaast.fasta

# <<<< adjust these two values to match the number of cores, GB of ram on the machine >>>>
NUM_CORES=12
MAX_MEM=50

MAX_BYTE_MEM=$(( $MAX_MEM*1024*1024*1024 ))
RUN_JAVA="/usr/local/java/jdk1.6.0_22_x64/bin/java -Xmx$MAX_BYTE_MEM"

$RUN_JAVA >/dev/null 2>/dev/null &
waitForJobs
# silently test if we have enough memory... but fail early if we don't

# <<<< these only need to be run once >>>>
python buildVAASTreference.py --in $REF_FASTA --out $REF_FASTA2 &
#python index_kgp.py --data $KGP_DATA_DIR --populations $KGP_POP_DIR &
MYDIR=`pwd`
#cd $SNPEFF_DIR
#$RUN_JAVA -jar snpEff.jar download -v hg19 &
#cd $MYDIR
# You also need to download dbNSFP2.0***.txt.gz and decompress it in SNPEFF_DIR,
# but you'll need to do that via web browser
waitForJobs

# <<<< reset this for a different project >>>>
TARGET_DIR=/raid1/alex/sequencing/cll/runs/exome_22_Oct_2012
TEMP=(`ls $TARGET_DIR/calls/all.*.finished.vcf`)
TARGET_VCFS=""
for s in ${TEMP[*]}
do
	v=${s##*/}
	TARGET_VCFS="$TARGET_VCFS ${v%.vcf}"
done

# I use this if statement to comment stuff out
if [[ 0 == 1 ]]
then
# ******** Actual jobs start here ********
rm -rf $TARGET_DIR/annotation
mkdir $TARGET_DIR/annotation

echo "Running snpEff"
echo "--------------"

rm -rf $TARGET_DIR/annotation/snpeff
mkdir $TARGET_DIR/annotation/snpeff
mkdir $TARGET_DIR/annotation/snpeff/logs

echo "Straight snpEff"
for s in ${TARGET_VCFS[*]}
do
	$RUN_JAVA -jar $SNPEFF_DIR/snpEff.jar \
		-c $SNPEFF_DIR/snpEff.config \
		-s $TARGET_DIR/annotation/snpeff/$s.summary.html \
		hg19 \
		$TARGET_DIR/calls/$s.vcf \
		>$TARGET_DIR/annotation/snpeff/$s.predbnsfp.vcf \
		2>$TARGET_DIR/annotation/snpeff/logs/$s.predbnsfp.err.log &
done
waitForJobs

echo "dbNSFP"
for s in ${TARGET_VCFS[*]}
do
	$RUN_JAVA -jar $SNPEFF_DIR/SnpSift.jar \
		dbnsfp \
		-a \
		$DBNSFP \
		$TARGET_DIR/annotation/snpeff/$s.predbnsfp.vcf \
		>$TARGET_DIR/annotation/snpeff/$s.pregwascat.vcf \
		2>$TARGET_DIR/annotation/snpeff/logs/$s.pregwascat.err.log &
done
waitForJobs

echo "GWAS catalog"
for s in ${TARGET_VCFS[*]}
do
	$RUN_JAVA -jar $SNPEFF_DIR/SnpSift.jar \
		gwasCat \
		$GWAS_CAT \
		$TARGET_DIR/annotation/snpeff/$s.pregwascat.vcf \
		>$TARGET_DIR/annotation/snpeff/$s.annotated.vcf \
		2>$TARGET_DIR/annotation/snpeff/logs/$s.annotated.err.log &
done
waitForJobs

echo "Running VAAST"
echo "-------------"

rm -rf $TARGET_DIR/annotation/vaast
mkdir $TARGET_DIR/annotation/vaast
mkdir $TARGET_DIR/annotation/vaast/logs

echo "Removing filtered variants"
for s in ${TARGET_VCFS[*]}
do
	$RUN_JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T SelectVariants \
		--variant $TARGET_DIR/calls/$s.vcf \
		-o $TARGET_DIR/annotation/vaast/$s.vcf \
		-R $REF_FASTA \
		--excludeFiltered \
		>$TARGET_DIR/annotation/vaast/logs/$s.selectVariants.log \
		2>$TARGET_DIR/annotation/vaast/logs/$s.selectVariants.err.log &
done
waitForJobs

fi

echo "Picking only the chromosomes VAAST likes (yuck... this is a problem!)"
for s in ${TARGET_VCFS[*]}
do
	python simplifyChromosomes.py \
		--in $TARGET_DIR/annotation/vaast/$s.vcf \
		--out $TARGET_DIR/annotation/vaast/$s.simple.vcf \
		>$TARGET_DIR/annotation/vaast/logs/$s.simplifyChromosomes.log \
		2>$TARGET_DIR/annotation/vaast/logs/$s.simplifyChromosomes.err.log &
done
waitForJobs

echo "Converting"
for s in ${TARGET_VCFS[*]}
do
	rm -rf $TARGET_DIR/annotation/vaast/$s
	mkdir $TARGET_DIR/annotation/vaast/$s
	mkdir $TARGET_DIR/annotation/vaast/$s/pre
	$VAAST_CONVERTER \
		--build hg19 \
		--format VCF \
		--path $TARGET_DIR/annotation/vaast/$s/pre/ \
		$TARGET_DIR/annotation/vaast/$s.simple.vcf \
		>$TARGET_DIR/annotation/vaast/logs/$s.vaast_converter.log \
		2>$TARGET_DIR/annotation/vaast/logs/$s.vaast_converter.err.log &
done
waitForJobs

echo "Sorting"
for s in ${TARGET_VCFS[*]}
do
	for i in $TARGET_DIR/annotation/vaast/$s/pre/*.gvf
	do
		j=${i##*/}
		$VAAST_SORT_GFF \
			--in_place \
			--perl_sort \
			--no_backup \
			$i \
			>$TARGET_DIR/annotation/vaast/logs/$s.${j%.gvf}.pre.vaast_sort_gff.log \
			2>$TARGET_DIR/annotation/vaast/logs/$s.${j%.gvf}.pre.vaast_sort_gff.err.log &
	done
done
waitForJobs

echo "Running VAT"
for s in ${TARGET_VCFS[*]}
do
	rm -rf $TARGET_DIR/annotation/vaast/$s/post
	mkdir $TARGET_DIR/annotation/vaast/$s/post
	for i in $TARGET_DIR/annotation/vaast/$s/pre/*.gvf
	do
		j=${i##*/}
		$VAT \
			--build hg19 \
			--fasta $REF_FASTA2 \
			--features $REF_GENE \
			$i \
			>$TARGET_DIR/annotation/vaast/$s/post/$j \
			2>$TARGET_DIR/annotation/vaast/logs/$s.${j%.gvf}.VAT.log &
	done
done
waitForJobs

echo "Sorting"
for s in ${TARGET_VCFS[*]}
do
	for i in $TARGET_DIR/annotation/vaast/$s/post/*.gvf
	do
		j=${i##*/}
		$VAAST_SORT_GFF \
			--in_place \
			--perl_sort \
			--no_backup \
			$i \
			>$TARGET_DIR/annotation/vaast/logs/$s.${j%.gvf}.post.vaast_sort_gff.log \
			2>$TARGET_DIR/annotation/vaast/logs/$s.${j%.gvf}.post.vaast_sort_gff.err.log &
	done
done
waitForJobs

echo "Running VST"
for s in ${TARGET_VCFS[*]}
do
	$VST \
		--ops 'U(0..5)' \
		$TARGET_DIR/annotation/vaast/$s/post/*.gvf \
		>$TARGET_DIR/annotation/vaast/$s.cdr &
done
waitForJobs

echo "Checking .cdr Quality"
for s in ${TARGET_VCFS[*]}
do
	$QUALITY_CHECK_PL \
		-sim 100000 \
		$VAAST_BACKGROUND \
		$TARGET_DIR/annotation/vaast/$s.cdr \
		>$TARGET_DIR/annotation/vaast/logs/$s.quality-check.log \
		2>$TARGET_DIR/annotation/vaast/logs/$s.quality-check.err.log &
done
waitForJobs

echo "Running VAAST"
for s in ${TARGET_VCFS[*]}
do
	$VAAST \
		--mode lrt \
		--outfile $TARGET_DIR/annotation/vaast/$s \
		--indel \
		$REF_GENE \
		$VAAST_BACKGROUND \
		$TARGET_DIR/annotation/vaast/$s.cdr \
		>$TARGET_DIR/annotation/vaast/logs/$s.VAAST.log \
		2>$TARGET_DIR/annotation/vaast/logs/$s.VAAST.err.log &
done
waitForJobs

if [[ 0 == 1 ]]
then
echo "Running calcStats.py"
echo "--------------------"

rm -rf $TARGET_DIR/annotation/calcStats
mkdir $TARGET_DIR/annotation/calcStats

for s in ${TARGET_VCFS[*]}
do
	python calcStats.py \
		--data $KGP_DATA_DIR \
		--in $TARGET_DIR/calls/$s.vcf \
		--out $TARGET_DIR/annotation/calcStats/$s.csv \
		>$TARGET_DIR/annotation/calcStats/$s.calcStats.log \
		2>$TARGET_DIR/annotation/calcStats/$s.calcStats.err.log &
		waitForJobs
done
fi
