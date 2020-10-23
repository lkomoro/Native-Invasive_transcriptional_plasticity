#!/bin/bash
#run transdecoder part 1 (predict longest ORFs)

#BSUB -q long
#BSUB -W 10:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e transdecoderORF_DS.err
#BSUB -oo transdecoderORF_DS.log

module load  transdecoder/5.5.0
TransDecoder.LongOrfs -t Delta_final_annotation_filtered_contigs_simple_headers.fasta
