#!/bin/bash

for file in /home/lkomoroske/DDIGRNASeqdata/test.folder/*.sam

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"` 

echo $sample

/usr/local/samtools-1.2/bin/samtools view -bh -@6 -F2048 "$sample"_mem.sam | \
/usr/local/samtools-1.2/bin/samtools sort -m 16G -@6 -O bam \
-T temporarysort -o "$sample"_PEsort.bam - \
> "$sample"_PE_sort.stdout 2> "$sample"_PE_sort.stderr

done