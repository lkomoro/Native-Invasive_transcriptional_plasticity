#!/bin/bash

#removing secondary reads
mkdir ../rm_2reads
for file in *.bam

do

echo $file

sample=`echo $file | cut -f6 -d "/" |cut -f1 -d "_" `

echo $sample

/usr/local/samtools-1.3/bin/samtools view -bh -@6 -F256 "$sample"_unfiltsort.bam | \
/usr/local/samtools-1.3/bin/samtools sort -m 16G -@6 -O bam \
-T temporarysort -o ../rm_2reads/"$sample"_PEsortunfilt.bam  \
> ../rm_2reads/"$sample".stdout 2> ../rm_2reads/"$sample".stderr

done
