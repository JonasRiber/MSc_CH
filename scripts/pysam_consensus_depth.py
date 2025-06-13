# pysam

# identify overlapping reads
# return: read.name o.start o.end

# call multinucleotide variants


import sys
import pysam
import csv


bam_file_path = sys.argv[1]
called_mismatches = sys.argv[2]
output_cb = sys.argv[3]


with pysam.AlignmentFile(bam_file_path, "rb") as bamfile, open(called_mismatches, 'r') as fcalls, open(output_cb, 'wb') as fcb:
	#header = '#CHROM\tSTART\tEND\tREF\tALT\tMTYPE\tN_CONSENSUS\n'.encode('utf-8')
	#fmnv.write(header)
	#
	nrows = 0
	reader = csv.reader(fcalls, delimiter='\t')
	next(reader)
	for row in reader:
		chromi, starti, endi, refi = str(row[0]), int(row[1]), int(row[2]), str(row[3])
		#if nrows > 100:
		#	break
		# consensus bases
		cb = 0
		# store paired read names in set
		read_set = set()
		# loop through all reads overlapping the position

		for read in bamfile.fetch(reference=chromi, start=starti, end=endi):
			
			fragment_length = abs(read.template_length)
			# initial quality control
			if(read.mapping_quality<60
				or read.is_unmapped
				or read.is_secondary
				or read.is_supplementary
				or not read.is_proper_pair
				or fragment_length>=300):
				continue
			# extract cigar statistics
			cigar_stats = read.get_cigar_stats()[0]
			# NM tag: number of mismatches per read
			NM = cigar_stats[10]
			
			if(NM > 10): # exclusion criteria
				continue
			# look for indels in cigar
			has_indel = sum(cigar_stats[1:4]) != 0 # I, D, N (BAM_CINS, BAM_CDEL, BAM_CREF_SKIP)
			if(has_indel):
				continue
			# look for clipping in cigar
			has_clip  = sum(cigar_stats[4:6]) != 0 # S, H (BAM_CSOFT_CLIP, BAM_CHARD_CLIP)
			if(has_clip):
				continue
			# store info on current read
			start = read.reference_start	# read start
			end   = read.reference_end		# read end
			name  = read.query_name			# read name
			# quality control
			if (start is None or end is None):
				continue

			positions = list(range(start, end)) # list positions in read
			# aligned sequence (i.e. not soft-clipped sequence)
			# This is a substring of query_sequence that excludes flanking bases that were soft clipped
			seq    = str(read.query_alignment_sequence).upper()
			refseq = str(read.get_reference_sequence()).upper()
			qual   = read.query_qualities


			for i in range(len(positions)):
				pos = positions[i]
				if pos == starti:
					if(qual[i]!=37):
						continue

					read_base = f"{name}_{pos}"
					if not read_base in read_set:
						read_set.add(read_base)
						cb += 1
					else:
						continue 	
			
			# for i in range(len(positions)):
			# 	pos = positions[i]
			# 	if pos == starti:
			# 		if(qual[i]!=37):
			# 			continue
					
			# 		cb += 1 # count base in position

			# 		base = seq[i]
			# 		read_base = f"{name}_{pos}_{base}"
			# 		print(read_base)
			# 		if not read_base in read_set:
			# 			read_set.add(read_base)
			# 		else:
			# 			cb -= 1 	#consensus so avoid counting multiple times
			# 		# 	cb += 1		# count when seen on second read in pair (consensus)
		nrows += 1
		record = f"{chromi}\t{starti}\t{endi}\t{refi}\t{cb}\n"
		#print(record)
		fcb.write(record.encode('utf-8'))

print("DONE")