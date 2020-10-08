#!/bin/bash
#scythe is simple adaptor trimmer

#module load scythe/c128b19 this is from monica, don't think we need it

#for details see github README file, also type scythe into command line and options will come up
#-q tells it it's quality encoding (default is Sanger, which is what ours is with Casava 1.8; earlier versions I think uses 'illumina', see Scythe documentation
#-n is a minimum match length argument; default is 5
#-M, --min-keep	filter sequnces less than or equal to this length (default: 35)
#-o, --output-file	output trimmed sequences file (default: stdout); your input file comes after this
#-p, --prior		prior (default: 0.300)
#-a path for your adaptor file

for file in *_goodreads.fastq

do

echo $file

sample=`echo $file | cut -f1 -d "_"`

echo $sample

scythe -M 30 \
-a /home/lkomoroske/DDIGRNASeqdata/NEB.adaptor.fasta \
-q sanger -n 0 \
-o "$sample"_R1_at.fastq \
"$sample"_combined_R1_goodreads.fastq \
> "$sample"_R1_scythe.stdout \
2> "$sample"_R1_scythe.stats.txt &

scythe -M 30 \
-a /home/lkomoroske/DDIGRNASeqdata/NEB.adaptor.fasta \
-q sanger -n 0 \
-o "$sample"_R2_at.fastq \
"$sample"_combined_R2_goodreads.fastq \
> "$sample"_R2_scythe.stdout \
2> "$sample"_R2_scythe.stats.txt

done

