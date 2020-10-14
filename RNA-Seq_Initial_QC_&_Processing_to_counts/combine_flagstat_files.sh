#!/bin/bash

for filename in *.flagstat; do
    echo "$filename"
    cat "$filename"
done > combined.flagstat.txt