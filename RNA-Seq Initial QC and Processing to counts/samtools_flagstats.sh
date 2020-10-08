#!/bin/bash


for file in /home/lkomoroske/DDIGRNASeqdata/All.BAM-FILES/*.bam

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"` 

echo $sample

/usr/local/samtools-1.3/bin/samtools flagstat "$sample"_PEsort.bam \
> "$sample"_PEsort.bam.flagstat 2> "$sample"_PEsort.bam.bai.stderr \

/usr/local/samtools-1.3/bin/samtools flagstat "$sample"_singleton_sort.bam \
> "$sample"_singletonsort.bam.flagstat 2> "$sample"_singletonsort.bam.stderr \

done


