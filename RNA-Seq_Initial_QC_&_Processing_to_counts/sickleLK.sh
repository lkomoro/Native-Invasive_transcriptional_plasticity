#!/bin/bash
#sickle is a windowed adaptive trimmer tool

#module load sickle/7667f147e6
#call for paired end with the 'sickle pe'
#then tell it how to find the forward (-f) and reverse (-r) of the pairs
# -t is the type (Casava 1.8 or higher is sanger)


for file in *R2_at.fastq

do

echo $file

sample=`echo $file | cut -f1 -d "_"`

echo $sample

sickle pe \
-f "$sample"_R1_at.fastq \
-r "$sample"_R2_at.fastq \
-t sanger -l 30 \
-o ../Mb-qt-files/"$sample"_R1_qt.fastq \
-p ../Mb-qt-files/"$sample"_R2_qt.fastq \
-s ../Mb-qt-files/"$sample"_single_qt.fastq \
> ../Mb-qt-files/"$sample"_sickle.stdout \
2> ../Mb-qt-files/"$sample"_sickle.stderr

done
