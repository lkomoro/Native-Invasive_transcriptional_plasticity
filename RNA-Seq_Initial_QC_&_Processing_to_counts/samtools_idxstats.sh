#!/bin/bash


for file in *.bam

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"`

echo $sample

/usr/local/samtools-1.3/bin/samtools idxstats "$sample"_PEsortv2.bam \
> "$sample"_PE.stdout 2> "$sample"_PE.stderr \

/usr/local/samtools-1.3/bin/samtools idxstats "$sample"_singleton_sortv2.bam \
> "$sample"_SE.stdout 2> "$sample"_SE.stderr \

done
