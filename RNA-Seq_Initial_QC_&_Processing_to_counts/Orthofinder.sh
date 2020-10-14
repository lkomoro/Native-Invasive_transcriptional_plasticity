#!/bin/bash
#run Orthofinder on translated fasta files


#BSUB -q long
#BSUB -W 72:00
#BSUB -R rusage[mem=16000]
#BSUB -n 12
#BSUB -R span[hosts=1]
#BSUB -e Orthofinder_full.err
#BSUB -oo Orthofinder_full.log

#run from: /project/uma_lisa_komoroske/Lisa/DDIG/Orthofinder_analyses/
    module load orthofinder/2.3.3
    singularity exec $ORTHOFINDERIMG orthofinder -f Orthofinder_2020

# SIMPLE USAGE:
# Run full OrthoFinder analysis on FASTA format proteomes in <dir>
#   orthofinder [options] -f <dir>
# 
# Add new species in <dir1> to previous run in <dir2> and run new analysis
#   orthofinder [options] -f <dir1> -b <dir2>
# 
# OPTIONS:
#  -t <int>          Number of parallel sequence search threads [Default = 40]
#  -a <int>          Number of parallel analysis threads [Default = 1]
#  -M <txt>          Method for gene tree inference. Options 'dendroblast' & 'msa'
#                    [Default = dendroblast]
#  -S <txt>          Sequence search program [Default = diamond]
#                    Options: blast
#  -A <txt>          MSA program, requires '-M msa' [Default = mafft]
#                    Options: mafft
#  -T <txt>          Tree inference method, requires '-M msa' [Default = fasttree]
#                    Options: fasttree
#  -s <file>         User-specified rooted species tree
#  -I <int>          MCL inflation parameter [Default = 1.5]
#  -x <file>         Info for outputting results in OrthoXML format
#  -p <dir>          Write the temporary pickle files to <dir>
#  -1                Only perform one-way sequence search
#  -X                Don't add species names to sequence IDs
#  -n <txt>          Name to append to the results directory
#  -o <txt>          Non-default results directory
#  -h                Print this help text
# 
# WORKFLOW STOPPING OPTIONS:
#  -op               Stop after preparing input files for BLAST
#  -og               Stop after inferring orthogroups
#  -os               Stop after writing sequence files for orthogroups
#                    (requires '-M msa')
#  -oa               Stop after inferring alignments for orthogroups
#                    (requires '-M msa')
#  -ot               Stop after inferring gene trees for orthogroups
# 
# WORKFLOW RESTART COMMANDS:
#  -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>
#  -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>
#  -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>

