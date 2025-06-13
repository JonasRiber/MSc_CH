#!/bin/bash


########################
### Merge call files ###
########################

### Author: Jonas Riber Jørgensen
### Date: 15-04-2025
### Description: exclude recurrent germline calls 


# input args
INPUT=$1
PATIENTID=$2
OUTPUT=$3

#Calls to be excluded
EXCLUDE=/path/to/gustavs/germline_calls/"$PATIENTID"/"$PATIENTID"_snpcalls_QUAL20_stack_geQbinom.tsv

#log file
log_file=data/"$PATIENTID"/"$PATIENTID"_filter_log.txt
before_filter=$(wc -l < "$INPUT")

mkdir -p tmp
SORTED_INPUT="tmp/${PATIENTID}_input_sorted.bed"
SORTED_EXCLUDE="tmp/${PATIENTID}_exclude_sorted.bed"

# Sort input and blacklist files
tail -n +2 "$INPUT" | sort -k1,1 -k2,2n > "$SORTED_INPUT"
tail -n +2 "$EXCLUDE" | sort -k1,1 -k2,2n > "$SORTED_EXCLUDE"


# Subtract exclusion regions	
bedtools subtract -a "$SORTED_INPUT" -b "$SORTED_EXCLUDE" > "$OUTPUT"


after_filter=$(wc -l < "$OUTPUT")
reduction=$(echo "scale=2; ($before_filter - $after_filter) / $before_filter * 100" | bc)
#write current entry
echo -e "$(date)\t| germline blacklisted calls |\t$before_filter\t$after_filter\t$reduction" >> "$log_file"

#clean up temp files
rm "$SORTED_INPUT"
rm "$SORTED_EXCLUDE"
