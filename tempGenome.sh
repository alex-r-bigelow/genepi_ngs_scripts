#!/bin/bash
# PHASE_START		"setup", "align", "build_bam", "sort_bam", "post_process", "call", "filter", "annotate"
#					Phase to begin with - most of the time it will be align unless you've added a new reference genome
#					you're restarting after a crash
export PHASE_START="realign"
# PHASE_STOP		"align", "build_bam", "sort_bam", "post_process", "call", "filter", "annotate", ""
#					Phase before which to end... "" will run the whole pipeline. Useful for debugging
export PHASE_STOP=""

# DATA_DIR			Directory. Should contain directories for each sample (named appropriately - these
#					names will be reused all the way past the call phase). Each subdirectory should
#					contain paired *.R1.*fastq.gz and *.R2.*fastq.gz (each lane should have exactly
#					the same name except for the R1 and R2)
export DATA_DIR=/raid1/alex/sequencing/cll/data/genomes

# TARGET_DIR		The directory in which to store files, run everything. You can safely experiment
#					with MIN_PRUNING and/or FILTER_MODE on the same directory; otherwise you should only
#					restart on existing data if something failed
CURRENT_DATE=`date +"%d_%b_%Y"`
# <<<< comment this to start a new run >>>>
CURRENT_DATE=14_Nov_2012
export TARGET_DIR=/raid1/alex/sequencing/cll/runs/genome_$CURRENT_DATE

# MIN_PRUNING		blank or an integer. If blank, unifiedGenotyper will be used for variant calling,
#					otherwise haplotypeCaller will be used with --minPruning MIN_PRUNING
#					(a low value for MIN_PRUNING will be very slow)
#export MIN_PRUNING=

# EXOME_OR_GENOME	"exome" or "genome"
#					will follow the appropriate parameters, depending on whether we're dealing with
#					exome/targeted or whole genome.
export EXOME_OR_GENOME="genome"

# EXOME_TARGETS		path to .bed file
#					If you're doing exome/targeted sequencing, you must provide a .bed file for
#					the targeted regions
# (set below)

# NUM_CORES			Number of processors on the machine running jobs
export NUM_CORES=12

# MAX_MEM			max GB RAM to allow these jobs to consume
export MAX_MEM=64

# FILTER_MODE		"ALL", "PURGED", or "PASS" - "ALL" will run annotations/summaries on all
#					variants, "PURGED" will run annotations on only on-target variants for
#					exome/targeted studies, and "PASS" will run annotations only on variants that
#					PASS all filters (except the DP Filter, which is applied by hand afterward)
export FILTER_MODE="PASS"

# VST_SET_OP		String to give VST for combining individuals
export VST_SET_OP="U(0..5)"

# Paths to the main directory of each app:
# JAVA
export JAVA=/usr/local/java/jdk1.6.0_22_x64

# BWA_DIR
export BWA_DIR=/raid1/sequencing/apps/alignment/bwa-0.6.2

# SAM_DIR
export SAM_DIR=/raid1/sequencing/apps/post_processing/samtools-0.1.18

# PICARD_DIR
export PICARD_DIR=/raid1/sequencing/apps/post_processing/picard-tools-1.78

# GATK_DIR
export GATK_DIR=/raid1/sequencing/apps/post_processing/GenomeAnalysisTK-2.3-4-g57ea19f

# SNPEFF_DIR
export SNPEFF_DIR=/raid1/sequencing/apps/annotation/snpEff_3_1

# VAAST
export VAAST=/raid1/sequencing/apps/analysis/VAAST/VAAST_RC_1.0.4/bin/VAAST

# VAAST_vaast_converter
export VAAST_vaast_converter=/raid1/sequencing/apps/analysis/VAAST/VAAST_RC_1.0.4/bin/vaast_tools/vaast_converter

# VAAST_VST
export VAAST_VST=/raid1/sequencing/apps/analysis/VAAST/VAAST_RC_1.0.4/bin/VST

# VAAST_vaast_sort_gff
export VAAST_vaast_sort_gff=/raid1/sequencing/apps/analysis/VAAST/VAAST_emailed/vaast_sort_gff

# VAAST_VAT
export VAAST_VAT=/raid1/sequencing/apps/analysis/VAAST/VAAST_RC_1.0.4/bin/VAT

# VAAST_QC
export VAAST_QC=/raid1/sequencing/apps/analysis/VAAST/VAAST_RC_1.0.4/bin/vaast_tools/quality-check.pl

# Paths to reference files:
REF_DIR=/raid1/sequencing/reference/gatk_bundle/gatk_bundle_04_Oct_2012
# --- Reference genomes ----
# REF_FASTA
export REF_FASTA=$REF_DIR/human_g1k_v37.fasta
# VAAST_FASTA		VAAST doesn't cooperate well with reference genomes other than what they give you;
#					I recommend making your own VAAST-compatible reference with buildVAASTreference.py,
#					but here's where you'd stick one they give you if you want to do that instead
# export VAAST_FASTA=/raid1/sequencing/reference/vaast/hg19/chrAll.fa
export VAAST_FASTA=$REF_DIR/human_g1k_v37.vaast.fasta

# --- GATK bundle ---
# REF_DBSNP_129
export REF_DBSNP_129=$REF_DIR/dbsnp_135.b37.excluding_sites_after_129.vcf

# REF_DBSNP
export REF_DBSNP=$REF_DIR/dbsnp_135.b37.vcf

# REF_MILLS
export REF_MILLS=$REF_DIR/Mills_and_1000G_gold_standard.indels.b37.sites.vcf

# REF_KGP
export REF_KGP=$REF_DIR/1000G_phase1.indels.b37.vcf

# REF_HAPMAP
export REF_HAPMAP=$REF_DIR/hapmap_3.3.b37.vcf

# REF_OMNI
export REF_OMNI=$REF_DIR/1000G_omni2.5.b37.sites.vcf

# (EXOME_TARGETS - see above)
export EXOME_TARGETS=$DATA_DIR/SeqCap_EZ_Exome_v3_b37_capture+300.bed

# ---- snpEff ----
# GWAS_CAT
export GWAS_CAT=/raid1/sequencing/reference/other/gwascatalog.txt
# DBNSFP
export DBNSFP=/raid1/sequencing/reference/other/dbNSFP2.0b3.txt

# ---- VAAST ----
# VAAST_BACKGROUND
export VAAST_BACKGROUND=/raid1/sequencing/reference/background/vaast/1KG_refGene_Dec2011_CGDiv_NHLBI_NoCall.cdr
# VAAST_FEATURES
export VAAST_FEATURES=/raid1/sequencing/reference/vaast/refGene_hg19.gff3

s=`basename $0`
/raid1/sequencing/apps/wrapper/genepi_ngs_scripts/best_practice_v4.sh >$s.log 2>$s.err.log