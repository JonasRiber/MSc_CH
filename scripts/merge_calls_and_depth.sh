#!/bin/bash

##################################
### Merge call and depth files ###
##################################

### Author: Jonas Riber Jørgensen
### Date: 25-02-2025
### Description: Merge mismtaches and depth from seperate files into one


#input args
INPUT_CALLS=$1
INPUT_DEPTH=$2

OUTPUT=$3


echo -e "#CHROM\tSTART\tEND\tREF\tALT\tMTYPE\tN_CON\tN_ALT\tN_ALL\tVAF" > "$OUTPUT"
awk '
BEGIN { FS = OFS = "\t"}                    # ensures tab seperation
FNR==NR {
    key = $1 FS $2 FS $3 FS $4                              #create a unique key using chr, start, end, ref
    data[key]=$5                                          #store the extra column from file2
    next
}
FNR > 1 {                                                   #skip first row (header)
    key = $1 FS $2 FS $3 FS $4
    depth = data[key]
    VAF = ($8 && depth) ? $8/depth : "."                     #calc VAF if both exists
    print $0, (key in data ? depth : "."), VAF      #append value from file2 if found, else add "."
}' "$INPUT_DEPTH" "$INPUT_CALLS" >> "$OUTPUT"





# bedtools intersect -wa -wb -a "$INPUT_CALLS" -b "$INPUT_DEPTH" | \
# awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >> "$OUTPUT"



### result from bedtools intersect

#| \
# awk 'BEGIN{FS=OFS="\t"}{VAF=$7/$12; print $1,$2,$3,$4,$5,$6,$7,$12,VAF}' 
# 1       2               3               4       5       6       7       8       9               10              11      12      (13)
# chr22   10662384        10662385        C       T       1       1       chr22   10662384        10662385        C       99




