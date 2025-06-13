#!/bin/bash

# Check if there are arguments (input files)
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_file1> <input_file2> ... <input_fileN>"
    exit 1
fi

#output file
OUTPUT="$1"

# list of input files
INPUT="${@:2}"


#collect heaer
head -1 $(echo $INPUT | cut -d' ' -f1) > "$OUTPUT"

#concatenate all input files and process them
cat "$@" | \
grep -v "#" | \
sort -k1,1V -k2,3n -k4,5 | \
bedtools groupby -g 1,2,3,4,5,6 -c 7,8 -o sum,sum >> "$OUTPUT"


#remove multi-allelic positions (picks the most mismatched)
tmpfile=$(mktemp)

header=$(head -n 1 "$OUTPUT")

awk -F'\t' '
NR > 1 {
    key = $1 FS $2 FS $3 FS $4;
    if (!(key in best) || $8 > best[key]) {
        best[key] = $8;
        line[key] = $0;
    }
}
END {
    for (k in line) {
        print line[k];
    }
}' "$OUTPUT" > "$tmpfile"


#write the final output
echo "$header" > "$OUTPUT"
cat "$tmpfile" >> "$OUTPUT"

rm -f "$tmpfile"