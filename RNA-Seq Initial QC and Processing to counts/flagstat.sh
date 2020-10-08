#!/bin/bash
#this is just combining all the bam stats files for each sample into one text file so can ditch the rest

for filename in *flagstat.txt; do
    echo "$filename"
    cat "$filename"
done > All_Mbunfilt_sortfltr_combined.flagstat.txt
