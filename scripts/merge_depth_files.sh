#!/bin/bash


#########################
### Merge depth files ###
#########################

### Author: Jonas Riber Jørgensen
### Date: 24-02-2025
### Description: Merged depth files from the mnv caller such that the coverages are summed together


# check args given

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <output_file> <input_files...>"
    exit 1
fi


#define args
# output file
OUTPUT="$1"

# list of input files
INPUT="${@:2}"

#merge all files, sort, and sum N_CONSENSUS (col7)
# cat: grab the inputs
# sort -k1,1 -k2,2n:
#   Sorts the input by the first column (#CHROM) and second column (START) numerically
# bedtools groupby -g 1,2,3,4 -c 7 -o sum:
#   -g 1,2,3,4: Groups by the columns (#CHROM, START, END, MTYPE)
#   -c 5: Specifies the 5th column (depth) for aggregation
#   -o sum: Sums the values in the 5th column (depth) for each group

cat $INPUT | sort -k1,1 -k2,2n -k3,3n | bedtools groupby -g 1,2,3,4 -c 5 -o sum > "$OUTPUT"

echo "Merged file saved as $OUTPUT"
