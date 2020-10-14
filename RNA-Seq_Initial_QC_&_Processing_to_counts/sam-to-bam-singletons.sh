#!/bin/bash

for file in /home/lkomoroske/DDIGRNASeqdata/All-singleton-SAM-FILES/*.sam

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"` 

echo $sample

/usr/local/samtools-1.2/bin/samtools view -bh -@10 -F2048 "$sample"_singleton_mem.sam | \
/usr/local/samtools-1.2/bin/samtools sort -m 16G -@10 -O bam \
-T temporarysort -o "$sample"_singleton_sort.bam - \
> "$sample"_singleton_sort.stdout 2> "$sample"_singleton_sort.stderr

done
