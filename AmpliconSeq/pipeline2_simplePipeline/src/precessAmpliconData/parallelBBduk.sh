#!/bin/bash

# This script removes Illumina's sequencing adapters from the 3' ends.
# It process either single-end or paired-end data.
# At the configuration you should add the path for BBMAP, folder where you
# have stored the raw reads, a .txt file listing the input file name.
# You also should define if the data is single or paired in 'ispaired' parameter.
# In the 'suffix' parameter, you should add all expcept the sample name that 
# are present in the filename; for example, in our case we have D1-01_SE.fastq
# for a single-end file with 'D1-01' as the sample name, and 'SE' as the suffix.

# Stop script on error
set -ue

# Configuration
BBMAP=/Users/vitorpavinato/Softwares/bbmap
INPUTS=batch2
OUTPUTS=OUTDIR
ids=ids_se.txt
ispaired=false
suffix_se=SE
suffix_r1=R1
suffix_r2=R2

# Stop script on error
set -ue

# Let's make gnu parallel nagging go away.
echo 'will cite' | parallel --citation 2> /dev/null

# Check if the directory to save the cleaned reads DO NOT exists
if [ ! -d ${OUTPUTS} ] 
then
    mkdir -p ${OUTPUTS}
fi

if "$ispaired";
then
	# Run BBDuk adapters clean on PE raw reads.
    cat ${ids} | parallel ${BBMAP}/bbduk.sh -Xmx3g in1=${INPUTS}/{}_${suffix_r1}.fastq in2=${INPUTS}/{}_${suffix_r2}.fastq\
										           out1=${OUTPUTS}/{}_cleaned_1.fq out2=${OUTPUTS}/{}_cleaned_2.fq \
										           ecco=t ktrim=r k=25 mink=11 \
										           rcomp=t tbo=true minlen=35 stats=${OUTPUTS}/{}_stats.txt \
										           ref=${BBMAP}/resources/adapters.fa hdist=1 refstats=${OUTPUTS}/{}_refstats.txt; 
										           
										           
else
	# Run BBDuk adapters clean on SE raw reads.
    cat ${ids} | parallel ${BBMAP}/bbduk.sh -Xmx3g in=${INPUTS}/{}_${suffix_se}.fastq \
										           out=${OUTPUTS}/{}_cleaned.fq \
										           ktrim=r k=25 mink=11 \
										           rcomp=t minlen=35 stats=${OUTPUTS}/{}_stats.txt \
										           ref=${BBMAP}/resources/adapters.fa hdist=1 refstats=${OUTPUTS}/{}_refstats.txt; 
fi

wait
