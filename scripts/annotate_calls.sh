
#!/bin/bash


########################
###  Annotate calls  ###
########################

### Author: Jonas Riber Jørgensen
### Date: 28-02-2025
### Description: Annotate called positions with context from specified databases (refGeneWithVer)




INPUT=$1
OUTPUT=$2

# subtract 1 from the END col
# this way start and end index is the same when there is only one nucleotide
# When 2 nucleotides it should increase by 1 (according to ANNOVAR indexing)



#make ANNOVAR input format and switch to 1-based indexing (+1 to start)
awk 'NR>1 {OFS="\t"; $2=$2+1; print $0}' "$INPUT" > "$INPUT".avinput


#run annovar - calls the table_annovar function on the input
perl annovar/table_annovar.pl "$INPUT".avinput \
	annovar/humandb/ \
	-buildver hg38 \
	-out "$INPUT"_annotated \
	-protocol refGeneWithVer,gnomad_genome,clinvar_20140702,cosmic70 \
	-operation g,f,f,f \
	-remove \
	-nastring . \
	--otherinfo


head -n 1 "$INPUT"_annotated.hg38_multianno.txt

#format the output
awk -F'\t' 'BEGIN {OFS="\t"; print "Chr", "Start", "End", "Ref", "Alt", "Func", "Gene", "Exonic_func", "AA_change", "Gnomad_all", "Clinvar", "Cosmic", "Mtype", "Cons_count", "Alt_count", "Total_count", "VAF", "Somatic_threshold"} 
NR>1 {print $1, $2, $3, $4, $5, $6, $7, $9, $10, $11, $19, $20, $21, $22, $23, $24, $25, $26}' "$INPUT"_annotated.hg38_multianno.txt > "$OUTPUT"

rm "$INPUT".avinput 
rm "$INPUT"_annotated.hg38_multianno.txt
 
#paste "$INPUT_annotated.hg38_multianno.txt" extra_info.tsv > final_annotated_output.tsv

echo "DONE: Annotation complete"


### ANNOVAR output
#1		2		3			4	 5	     6                       7                       8                               9                              10                          11                    12         13                        14  15              
#Chr     Start   End     Ref     Alt     Func.refGeneWithVer     Gene.refGeneWithVer     GeneDetail.refGeneWithVer       ExonicFunc.refGeneWithVer      AAChange.refGeneWithVer  Otherinfo1      Otherinfo2      Otherinfo3      Otherinfo4      Otherinfo5
#1,2,3,4,5,6,7,9,11,12,13,14,15

#1		 2		 3		 4	 	 5	     6                       7                       8                               9                               10                      11                    	 12         		   13                      14  		               15                      16                      17                      18                      19                      20              21            22              23        
#Chr     Start   End     Ref     Alt     Func.refGeneWithVer     Gene.refGeneWithVer     GeneDetail.refGeneWithVer       ExonicFunc.refGeneWithVer       AAChange.refGeneWithVer gnomAD_genome_ALL       gnomAD_genome_AFR     gnomAD_genome_AMR       gnomAD_genome_ASJ       gnomAD_genome_EAS       gnomAD_genome_FIN       gnomAD_genome_NFE       gnomAD_genome_OTH       clinvar_20140702        cosmic70        Otherinfo1    Otherinfo2      Otherinfo3      Otherinfo4      Otherinfo5      Otherinfo6

