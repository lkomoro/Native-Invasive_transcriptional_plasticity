#!/bin/bash

#module load bwa/0.6.2
#module load samtools/1.2

cut -f1 -d " " ../Reference_fasta_files/M_beryllina_05_2013_61k.fasta > bwa_ref.fasta
#/nfs/data/export/mbritton/general_use/fasta1line.pl temp.fasta bwa_ref.fasta

bwa index -a is bwa_ref.fasta \
> bwa_index.stdout 2> bwa_index.stderr

#note I removed the -I option below adapting this script from BWA aln to BWA mem because BWA mem assumes phred+33
for file in /home/lkomoroske/DDIGRNASeqdata/Mb-qt-files/*R1_qt.fastq

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"` 

echo $sample


bwa mem -M -t 16 bwa_ref.fasta "$sample"_R1_qt.fastq "$sample"_R2_qt.fastq \
> "$sample"_mem.sam
2> "$sample"_PE_bwa_mem.stderr &


done
