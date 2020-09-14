## **Targeted Amplicon: Reference Sequences**

The file Belgica_antarctica_amplicons.fasta contains the refence sequence
of each targeted amplicon. It can be used as a reference sequence for 
read alignment. We used this sequences two in the pipeline: 1) to filter
out bad sequences within clustered sequences, 2) as a reference sequence
for read alignment.

## **Each amplicon sequence header contains the following values, separeted by "|"**
1. the amplicon name;
2. the scaffold name;
3. the interval [start-end] of the amplicon in the scaffold coordinates;
4. a list of SNPs (separated by ",").
