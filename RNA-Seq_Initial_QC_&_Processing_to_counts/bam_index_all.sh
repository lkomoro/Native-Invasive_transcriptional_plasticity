#!/bin/bash

for file in *.bam

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"` 

echo $sample

/usr/local/samtools-1.2/bin/samtools index "$sample"_PEsortv2.bam 
/usr/local/samtools-1.2/bin/samtools index  "$sample"_singleton_sortv2.bam 

done
