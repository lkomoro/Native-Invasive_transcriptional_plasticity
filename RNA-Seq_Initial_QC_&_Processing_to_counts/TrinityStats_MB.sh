#!/bin/bash
#run TrinityStats

#BSUB -q long
#BSUB -W 10:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e TrinityStats.err
#BSUB -oo TrinityStats.log

module load trinity/2.8.5
/share/pkg/trinity/2.8.5/util/TrinityStats.pl ./M_beryllina_05_2013_61k.fasta.transdecoder_dir/longest_orfs.cds
