#!/bin/bash
#run transdecoder Part 3

#BSUB -q long
#BSUB -W 10:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e transdecoderPredict_DS.err
#BSUB -oo transdecoderPredict_DS.log

module load  transdecoder/5.5.0
module load R/3.6.0
module load gcc/8.1.0
TransDecoder.Predict -t Delta_final_annotation_filtered_contigs_simple_headers.fasta
