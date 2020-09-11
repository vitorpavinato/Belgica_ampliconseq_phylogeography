#!/bin/bash

# This program performs read merging of raw or cleaned reads, insert size 
# calculation and adapter identification. It process only paired-end reads.
# At the configuration you should add the path for BBMAP, folder where you
# have stored the raw reads, a .txt file listing the input file name.
# In the 'suffix' parameter, you should add all expcept the sample name that 
# are present in the filename; for example, in our case we have D1-01_cleaned_1.fastq
# and D1-01_cleaned_2.fastq for first and second reads. The suffixes should be:
# 'cleaned_1' for the first reads and 'cleaned_2' for the second reads.

# Configuration
BBMAP=/Users/vitorpavinato/Softwares/bbmap
INPUTS=raw-cleaned
OUTPUTS=OUTDIR
ids=ids_pe.txt
suffix_r1=cleaned_1
suffix_r2=cleaned_2

# Stop script on error
set -ue

# Let's make gnu parallel nagging go away.
echo 'will cite' | parallel --citation 2> /dev/null

# Make a directory for the results.
mkdir -p ${OUTPUTS}

# Run BBmerge merge and adapters check on the raw reads.
cat ${ids} | parallel ${BBMAP}/bbmerge.sh -Xmx3g in1=${INPUTS}/{}_${suffix_r1}.fq \
                                           in2=${INPUTS}/{}_${suffix_r2}.fq \
                                           merge=t \
                                           out=${OUTPUTS}/{}_merged.fq \
                                           outu1=${OUTPUTS}/{}_unmerged_1.fq \
                                           outu2=${OUTPUTS}/{}_unmerged_2.fq \
                                           ihist=${OUTPUTS}/{}_ihist.txt \
                                           outa=${OUTPUTS}/{}_adapters.fa ; 


wait
 
