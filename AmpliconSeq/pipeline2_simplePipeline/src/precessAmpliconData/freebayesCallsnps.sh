#!/bin/bash

# This program performs the variant discovery on BAM files with Freebayes 

# Configuration
BAMLIST=bamlist.txt
mergebams=false
OUTPUTS_1=MergedBams
OUTPUTS_2=CalledVariants
INPUTS=sorted_merged.bam
REF=belgicaScaffold.fa

# Stop on errors. Print the commands.
set -uex

# Make a directory for the results.
# Check if the directory to save the cleaned reads DO NOT exists

if "$mergebams";
then
	if [ ! -d ${OUTPUTS_1} ] 
	then
		mkdir -p ${OUTPUTS_1}
	fi
	
	#Merge BAMs
	samtools merge -f -b ${BAMLIST} ${OUTPUTS_1}/sorted_merged.bam
	
	if [ ! -d ${OUTPUTS_2} ] 
	then
		mkdir -p ${OUTPUTS_2}
	fi
	
	# Run Freebayes
	freebayes --fasta-reference ${REF} ${OUTPUTS_1}/sorted_merged.bam > ${OUTPUTS_2}/raw_variants.vcf
	wait
else
	freebayes --fasta-reference ${REF} ${OUTPUTS_1}/sorted_merged.bam > ${OUTPUTS_2}/raw_variants.vcf
	wait

fi
