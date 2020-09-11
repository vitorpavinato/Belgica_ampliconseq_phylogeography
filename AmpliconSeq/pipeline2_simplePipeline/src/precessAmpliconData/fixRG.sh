#!/bin/bash

# This scripts fix the reading groups of read mapping files from BWA-MEM 
# and prepare the BAM files for subsequent analysis with GATK.

# Configuration
#PICARD=/usr/local/Cellar/picard-tools/2.23.3/libexec/
ids=ids_pe.txt
INPUTS=aligned_batch1
OUTPUTS=alignedRG_batch1

# Stop on errors. Print the commands.
set -uex

# Make a directory for the results.
# Check if the directory to save the cleaned reads DO NOT exists
if [ ! -d ${OUTPUTS} ] 
then
    mkdir -p ${OUTPUTS}
fi

# Run PICARD to add Reading Groups - it is a for loop that uses the information provided in ids.
for i in $(cat ${ids});
do
	#java -Xmx3g -jar ${PICARD}/picard.jar 
	picard SortSam VALIDATION_STRINGENCY=LENIENT I=${INPUTS}/${i}.bam O=${OUTPUTS}/${i}_sorted.bam SORT_ORDER=coordinate;
	wait

	#java -Xmx3g -jar ${PICARD}/picard.jar 
	picard AddOrReplaceReadGroups I=${OUTPUTS}/${i}_sorted.bam \
                                          RGID=${i} RGPU=None RGSM=${i} RGPL=Illumina RGLB=${i} O=${OUTPUTS}/${i}_sortedrg.bam
                                          
    samtools index ${OUTPUTS}/${i}_sortedrg.bam ${OUTPUTS}/${i}_sortedrg.bai

done
