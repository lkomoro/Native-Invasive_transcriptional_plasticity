#!/bin/bash

mkdir ../bamfiles
for file in /home/lkomoroske/DDIG/Mb_qt_files_pairs/*.sam

do

echo $file

sample=`echo $file | cut -f6 -d "/" |cut -f1 -d "_" `

echo $sample

/usr/local/samtools-1.3/bin/samtools view -bh -@6 -F2048 "$sample"_unfiltref.sam | \
/usr/local/samtools-1.3/bin/samtools sort -m 16G -@6 -O bam \
-T temporarysort -o "$sample"_unfiltsort.bam \
> "$sample"_sort.stdout 2> "$sample"_sort.stderr

done
