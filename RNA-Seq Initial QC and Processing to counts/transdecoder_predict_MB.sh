#!/bin/bash
#run transdecoder Part 3

#BSUB -q long
#BSUB -W 10:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e transdecoderPredict.err
#BSUB -oo transdecoderPredict.log

module load  transdecoder/5.5.0
module load R/3.6.0
module load gcc/8.1.0
TransDecoder.Predict -t M_beryllina_05_2013_61k.fasta --single_best_only

#  Transdecoder.LongOrfs|http://transdecoder.github.io> - Transcriptome Protein Prediction
#
#
#  Required:
#
#   -t <string>                            transcripts.fasta
#
#  Common options:
#
#
#   --retain_long_orfs_mode <string>        'dynamic' or 'strict' (default: dynamic)
#                                        In dynamic mode, sets range according to 1%FDR in random sequence of same GC content.
#
# 
#   --retain_long_orfs_length <int>         under 'strict' mode, retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence 
#                                         marks it as coding (default: 1000000) so essentially turned off by default.)
#
#   --retain_pfam_hits <string>            domain table output file from running hmmscan to search Pfam (see transdecoder.github.io for info)     
#                                        Any ORF with a pfam domain hit will be retained in the final output.
# 
#   --retain_blastp_hits <string>          blastp output in '-outfmt 6' format.
#                                        Any ORF with a blast match will be retained in the final output.
#
#   --single_best_only                     Retain only the single best orf per transcript (prioritized by homology then orf length)
#
#   --output_dir | -O  <string>            output directory from the TransDecoder.LongOrfs step (default: basename( -t val ) + ".transdecoder_dir")
#
#   -G <string>                            genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia, ...)
#
#   --no_refine_starts                     start refinement identifies potential start codons for 5' partial ORFs using a PWM, process on by default.
#
##  Advanced options
#
#    -T <int>                            Top longest ORFs to train Markov Model (hexamer stats) (default: 500)
#                                        Note, 10x this value are first selected for removing redundancies,
#                                        and then this -T value of longest ORFs are selected from the non-redundant set.
#  Genetic Codes
#
#
#   --genetic_code <string>                Universal (default)
#
#        Genetic Codes (derived from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
#
#
# Acetabularia
# Candida
# Ciliate
# Dasycladacean
# Euplotid
# Hexamita
# Mesodinium
# Mitochondrial-Ascidian
# Mitochondrial-Chlorophycean
# Mitochondrial-Echinoderm
# Mitochondrial-Flatworm
# Mitochondrial-Invertebrates
# Mitochondrial-Protozoan
# Mitochondrial-Pterobranchia
# Mitochondrial-Scenedesmus_obliquus
# Mitochondrial-Thraustochytrium
# Mitochondrial-Trematode
# Mitochondrial-Vertebrates
# Mitochondrial-Yeast
# Pachysolen_tannophilus
# Peritrich
# SR1_Gracilibacteria
# Tetrahymena
# Universal
#
#  --version                           show version (5.5.0)
#
#########################################################################################
