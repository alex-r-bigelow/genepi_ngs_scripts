EXOME_OR_GENOME=exome

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

DATA_DIR=/raid1/alex/sequencing/cll/data/$EXOME_OR_GENOME\s
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
EXOME_TARGETS=$DATA_DIR/SeqCap_EZ_Exome_v3_b37_capture.bed

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
CURRENT_DATE=22_Oct_2012
TARGET_DIR=/raid1/alex/sequencing/cll/runs/$EXOME_OR_GENOME\_$CURRENT_DATE

VCF_NAME=""
if [ -z $1 ]
then
	VCF_NAME="unifiedGenotyper"
else
	VCF_NAME="haplotypeCaller$1"
fi

waitForJobs