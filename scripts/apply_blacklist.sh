#!/bin/bash
# script to remove positions within blacklists and also remove initial low mismatched positions
# 


#check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <inputfile>"
    exit 1
fi

INPUT="$1"


# merge blacklists and subtract from variant calls
cat backup/blacklists/*.bed | cut -f1-3 | \
bedtools sort -i stdin | bedtools merge | \
bedtools subtract -a "$INPUT" -b stdin #|
#awk '$8 > 3'

