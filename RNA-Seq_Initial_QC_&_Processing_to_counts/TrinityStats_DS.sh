#!/bin/bash
#run TrinityStats

#BSUB -q long
#BSUB -W 10:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e TrinityStats_DS.err
#BSUB -oo TrinityStats_DS.log

module load trinity/2.8.5
/share/pkg/trinity/2.8.5/util/TrinityStats.pl ./Delta_final_annotation_filtered_contigs_simple_headers.fasta.transdecoder_dir/longest_orfs.cds
