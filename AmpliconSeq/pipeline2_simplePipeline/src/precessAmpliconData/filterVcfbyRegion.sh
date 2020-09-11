#!/bin/bash

# This scripts performs a SNP filtering based on an interval or on a list
# of SNPs (actually it is an interval of 1nt). 

# Configuration
VCFLIB=/Users/vitorpavinato/Softwares/vcflib/bin/
decompose=false
invert=true
INPUT=variants_outside_amplicon.vcf
OUTPUT_1=variants_outside_outsidePrimers.vcf
OUTPUT_2=NULL
INTVL=fluidigm_primers_interval_inGenome.bed

# Stop on errors. Print the commands.
set -uex

if "$decompose";
then
	${VCFLIB}/vcfallelicprimitives ${INPUT} > ${OUTPUT_1}
	
	if "$invert";
	then
		${VCFLIB}/vcfintersect -b ${INTVL} -v ${OUTPUT_1} > ${OUTPUT_2}
	
	else
		${VCFLIB}/vcfintersect -b ${INTVL} ${OUTPUT_1} > ${OUTPUT_2}

	fi # end of if

else

	if "$invert";
	then
		${VCFLIB}/vcfintersect -b ${INTVL} -v ${INPUT} > ${OUTPUT_1}
	
	else
		${VCFLIB}/vcfintersect -b ${INTVL} ${INPUT} > ${OUTPUT_1}

	fi # end of if

fi # end of if 

