#!/bin/bash
#this is just combining all the bam stats files for each sample into one text file that is tab delimited so one sample per row so can then easily manipulate to summarize mapping stats etc

for filename in *unfilt.flagstat.txt; do
    cat "$filename"|tr '\n' '\t'
    echo "$filename"
done > All_Mbunfilt_combined.flagstat_tabbed.txt
