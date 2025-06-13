#!/bin/sh

annotateQbinom() {
	if [[ $# -eq 0 ]] ; then
		echo "No arguments supplied. Two needed: 'p_threshold' and 'snpcalls' "
		exit 1
	elif [[ $# -eq 1 ]] ; then
		echo "Only 1 argument supplied. Two needed: 'p_threshold' and 'snpcalls' "
		exit 1
	elif [[ $# -ge 3 ]] ; then
		echo "Too many arugments supplied. Two needed: 'p_threshold' and 'snpcalls' "
	#elif ! [ -f $2 ] ; then
	#	echo "File does not exist: $2"
	fi

	# takes as input:
	# $1: p-value threshold
	# $2: snpcalls
	# define look-up table of qbinom thresholds
	p_threshold=$1
	calls=$2
	qbinom_table=/faststorage/project/c2i-colon/WorkSpaces/riber16/Master/mnv_workflow/gustav_scripts/qbinom_somatic_by_depth.tsv

	awk -v p_threshold=$p_threshold '
	# 
	BEGIN{ FS=OFS="\t"; # tab separated in- and outfile
	# create array of p-value thresholds (keys) and according column ids (values)
	P["0.1"]=2; P["0.05"]=3; P["0.01"]=4; P["1e-3"]=5;	P["1e-4"]=6;
	P["1e-5"]=7; P["1e-6"]=8; P["1e-7"]=9; P["1e-8"]=10; P["1e-9"]=11; P["1e-10"]=12}

	# build array of lower tail thresholds given the depth and p-value
	FNR>1 && FNR==NR {  				# skip header and operate on first input file
		Q[$1] = $P[p_threshold] 		# build array
	}

	FNR!=NR {
		print $0, Q[$9]
	}
	' $qbinom_table $calls
}

### Adjust the columns if used
filterQbinom() {
	awk '$5 <= $NF {print $0}' $1
}

filterINDELs() {
	awk 'length($3) == 1 && length($4) == 1 {print $0}' $1
}


cat $1 \
| annotateQbinom $2 - 
#| filterQbinom | filterINDELs

# how many single- or double-read supported calls or there?
#| awk '$5<=2{print $0}' #| wc -l

# compute mean min and max VAF
# | awk 'BEGIN{min=1;max=0} $7>max{max=$7} $7<min{min=$7} {vaf+=$7} END{print vaf/FNR, min, max}'


