#!/bin/bash

# This program performs the read mapping using BWA-MEM and prepare the 
# BAM files for subsequent analysis with GATK.

# Configuration
#PICARD=/usr/local/Cellar/picard-tools/2.20.8/libexec/
ids=ids_se.txt
ifjoin=false
INPUTS=raw_cleaned_batch2
OUTPUTS=OUTDIR
suffix1=cleaned
#suffix1=unmerged
#suffix2=merged
#REF=belgicaScaffold.fa

# Stop on errors. Print the commands.
set -uex

# Let's make gnu parallel nagging go away.
echo 'will cite' | parallel --citation 2> /dev/null

# Make a directory for the results.
# Check if the directory to save the cleaned reads DO NOT exists
if [ ! -d ${OUTPUTS} ] 
then
    mkdir -p ${OUTPUTS}
fi

# Make BWA reference index
#bwa index -p referenceidx ${REF}

if ${ifjoin};
then
   # Run BWA MEM to map 1st reads of unmerged PE reads to the reference.
   # Unmerged second PE reads were to bad to be used.
   cat ${ids} | time parallel -j+0 --eta 'bwa mem -aM -t 1 -R "@RG\tID:{}\tSM:{}\tPU:None\tPL:Illumina\tLB:{}" -o '${OUTPUTS}'/{}_'${suffix1}'.sam referenceidx '${INPUTS}'/{}_'${suffix1}'_1.fq'
   wait
   
   # Run SAMTOOLS convert SAM to BAM and sort BAM file.
   cat ${ids} | time parallel -j+0 --eta samtools view -bS ${OUTPUTS}/{}_${suffix1}.sam -o ${OUTPUTS}/{}_${suffix1}.bam
   wait
   
   # Run BWA MEM to map single reads to the reference.
   cat ${ids} | time parallel -j+0 --eta 'bwa mem -aM -t 1 -R "@RG\tID:{}\tSM:{}\tPU:None\tPL:Illumina\tLB:{}" -o '${OUTPUTS}'/{}_'${suffix2}'.sam referenceidx '${INPUTS}'/{}_'${suffix2}'.fq'
   wait
   
   # Run SAMTOOLS convert SAM to BAM and sort BAM file.
   cat ${ids} | time parallel -j+0 --eta samtools view -bS ${OUTPUTS}/{}_${suffix2}.sam -o ${OUTPUTS}/{}_${suffix2}.bam
   wait
   
   # Merge BAM files.
   cat ${ids} | time parallel -j+0 --eta samtools merge -f ${OUTPUTS}/{}.bam ${OUTPUTS}/{}_${suffix2}.bam ${OUTPUTS}/{}_${suffix1}.bam
   wait
   
   # Run Picard SortSam.
   #cat ${ids} | time parallel -j+0 --eta java -Xmx1g -jar ${PICARD}/picard.jar SortSam VALIDATION_STRINGENCY=LENIENT I=${OUTPUTS}/{}.bam O=${OUTPUTS}/{}_sorted.bam SORT_ORDER=coordinate
   #wait
   
   # Run Picard AddOrReplaceReadGroup.
   #cat ${ids} | time parallel -j+0 --eta java -Xmx1g -jar ${PICARD}/picard.jar AddOrReplaceReadGroups I=${OUTPUTS}/{}_sorted.bam \
                                                              RGID={} RGPU=None RGSM={} RGPL=Illumina RGLB={} O=${OUTPUTS}/{}_sortedrg.bam
   #wait
          
   # Run SAMTOOLS sort BAM.
   cat ${ids} | time parallel -j+0 --eta samtools sort -o ${OUTPUTS}/{}_sorted.bam ${OUTPUTS}/{}.bam 
   wait
   
   # Run SAMTOOLS index BAM.
   cat ${ids} | time parallel -j+0 --eta samtools index ${OUTPUTS}/{}_sorted.bam ${OUTPUTS}/{}_sorted.bai
   wait
   
   # Remove intermediate SAM files
   cat ${ids} | time parallel -j+0 --eta rm ${OUTPUTS}/{}_${suffix1}.sam
   wait
   
   cat ${ids} | time parallel -j+0 --eta rm ${OUTPUTS}/{}_${suffix2}.sam
   wait
   
   # Remove intermediate BAM files.
   cat ${ids} | time parallel -j+0 --eta rm ${OUTPUTS}/{}_${suffix1}.bam
   wait
   
   cat ${ids} | time parallel -j+0 --eta rm ${OUTPUTS}/{}_${suffix2}.bam
   wait
  
else
    
   # Run BWA MEM to map paired reads to the reference transcriptome.
   cat ${ids} | time parallel -j+0 --eta 'bwa mem -aM -t 1 -R "@RG\tID:{}\tSM:{}\tPU:None\tPL:Illumina\tLB:{}" -o '${OUTPUTS}'/{}_'${suffix1}'.sam referenceidx '${INPUTS}'/{}_'${suffix1}'.fq'
   wait
   
   # Run SAMTOOLS convert SAM to BAM and sort BAM file.
   cat ${ids} | time parallel -j+0 --eta samtools view -bS ${OUTPUTS}/{}_${suffix1}.sam -o ${OUTPUTS}/{}.bam
   wait
   
   # Run Picard SortSam.
   #cat ${ids} | time parallel -j+0 --eta java -Xmx1g -jar ${PICARD}/picard.jar SortSam VALIDATION_STRINGENCY=LENIENT I=${OUTPUTS}/{}.bam O=${OUTPUTS}/{}_sorted.bam SORT_ORDER=coordinate
   #wait
   
   # Run Picard AddOrReplaceReadGroup.
   #cat ${ids} | time parallel -j+0 --eta java -Xmx4g -jar ${PICARD}/picard.jar AddOrReplaceReadGroups I=${OUTPUTS}/{}_sorted.bam \
                                                               RGID={} RGPU=None RGSM={} RGPL=Illumina RGLB={} O=${OUTPUTS}/{}_sortedrg.bam
   #wait
          
   # Run SAMTOOLS sort BAM.
   cat ${ids} | time parallel -j+0 --eta samtools sort -o ${OUTPUTS}/{}_sorted.bam ${OUTPUTS}/{}.bam 
   wait
   
   # Run SAMTOOLS index BAM.
   cat ${ids} | time parallel -j+0 --eta samtools index ${OUTPUTS}/{}_sorted.bam ${OUTPUTS}/{}_sorted.bai
   wait
   
   # Remove intermediate SAM files.
   cat ${ids} | time parallel -j+0 --eta rm ${OUTPUTS}/{}_${suffix1}.sam
   wait
   
    
fi ## end of "ifjoin"

