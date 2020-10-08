#!/bin/bash
#run transdecoder part 1 (predict longest ORFs)

#BSUB -q long
#BSUB -W 10:00
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -e transdecoderORF.err
#BSUB -oo transdecoderORF.log

module load  transdecoder/5.5.0
TransDecoder.LongOrfs -t ../M_beryllina_05_2013_61k.fasta

########################################################################################
#
#  Transdecoder.LongOrfs|http://transdecoder.github.io> - Transcriptome Protein Prediction
#
#
#  Required:
#
#    -t <string>                            transcripts.fasta
#
#  Optional:
#
#   --gene_trans_map <string>              gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return> ) 
#
#   -m <int>                               minimum protein length (default: 100)
# 
#   -G <string>                            genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)
#  
#   -S                                     strand-specific (only analyzes top strand)
#
#   --output_dir | -O  <string>            path to intended output directory (default:  basename( -t val ) + ".transdecoder_dir")
#
#   --genetic_code <string>                Universal (default)
#
#        Genetic Codes (derived from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
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
#
#   --version                              show version tag (5.5.0)
#
#########################################################################################


