#!/bin/bash

cut -f1 -d " " ../Reference_fasta_files/Delta_final_annotation_filtered_contigs_simple_headers.fasta > bwa_ref.fasta

bwa index -a is bwa_ref.fasta \
> bwa_index.stdout 2> bwa_index.stderr

for file in /home/lkomoroske/DDIGRNASeqdata/issuessamples/*single_qt.fastq

do

echo $file

sample=`echo $file | cut -f1 -d "_" | cut -f6 -d "/"` 

echo $sample


bwa mem -M -t 16 bwa_ref.fasta "$sample"_single_qt.fastq \
> "$sample"_singleton_mem.sam
2> "$sample"_singleton_bwa_mem.stderr &


done
