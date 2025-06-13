# call multinucleotide variants

import sys
import pysam
import time

# from datetime import datetime
# from collections import Counter

start_time = time.time()

# args
bam_file_path = sys.argv[1]
output_mnv = sys.argv[2]
# input_chrom = sys.argv[3]
bed_path_file = sys.argv[3]

# small test:
# test - tmp/PT_ID_chr22_38640000_38650000.bam

# consensus script (pysam_mnv.py):
# No. of read pairs with consensus MNVs 311
# No. of MNVs (consensus): 53
# No. of mismatches (consensus) 559

# All script (this):
# No. of read pairs with consensus MNVs: 311
# No. of MNVs (total): 76
# No. of MNVs (consensus): 53
# No. of mismatches (total): 1034
# No. of mismatches (consensus) 559
# Consensus/total mismatch ratio: 0.5406189555125726
# Runtime: 0.0272 seconds


#note: it seems that the script collecting all mnvs is also 
# collecting the same consensus mnvs 


# larger test:
# chr2

# All script (this):
# No. of read pairs with consensus MNVs: 3201591
# No. of MNVs (total): 794689
# No. of MNVs (consensus): 505336
# No. of mismatches (total): 5770201
# No. of mismatches (consensus) 3582156
# Consensus/total mismatch ratio: 0.6208026375511009
# Saved consensus called MNVs to tmp/C06051_all_mm.tsv
# Script runtime: 324.9826 seconds


#parse bedfile
def parse_bed_file(bed_file_path):
    regions = []
    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            if line.strip() and not line.startswith('#'):
                fields = line.strip().split()
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                regions.append((chrom, start, end))
    return regions

regions = parse_bed_file(bed_path_file)



with pysam.AlignmentFile(bam_file_path, "rb") as bamfile, open(output_mnv, 'wb') as fmnv:
	header = '#CHROM\tSTART\tEND\tREF\tALT\tMTYPE\tN_CONSENSUS\tN_MISMATCHES\n'.encode('utf-8')
	fmnv.write(header)

	frag_dict = dict()	# fragment length counter
	read_c = 0			# read counter
	base_c = 0			# base counter
	read_set = set()	# store paired read names in set
	read_set_nonconsens = set()	# store paired read names in disagreement

	# store all MNVs with counts in dict
	mnv_dict_consensus = dict()
	mnv_dict = dict()
	mnv_set = set()

	mnv_dict_nonconsensus = dict()
	mnv_alt_dict = dict()
	mnv_no_alt_set = set()
	
	#loop over regions
	counter = 0
	nr_regions = len(regions)
	for chrom, start, end in regions:
		counter += 1
		print(f"Processing region {counter}/{nr_regions}: {chrom}:{start}-{end}")	
		
		# loop over reads in bamfile
		for read in bamfile.fetch(reference=chrom, start=start, end=end):

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
			NM = cigar_stats[10]	# NM tag: number of mismatches per read
			
			# if(NM != 1):
			if(NM == 0 or NM > 10): # exclusion criteria
				continue
			# look for indels in cigar
			has_indel = sum(cigar_stats[1:4]) != 0 	# I, D, N (BAM_CINS, BAM_CDEL, BAM_CREF_SKIP)
			if(has_indel):
				continue
			# look for clipping in cigar
			has_clip  = sum(cigar_stats[4:6]) != 0 	# S, H (BAM_CSOFT_CLIP, BAM_CHARD_CLIP)
			if(has_clip):
				continue
			
			# store info on current read
			chrom = read.reference_name		# chrom
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
			
			try:
				# extract position in read at mismatches
				mis_id = [i for i in range(len(seq)) if seq[i] != refseq[i]]
			except:
				# print(f'{name} problem')
				continue

			try:
				poss, refs, seqs = zip(*[
					(positions[i], refseq[i], seq[i])
					for i in mis_id
					if qual[i] == 37 and seq[i] != 'N' and refseq[i] != 'N'
				])
			except ValueError: # Handles case where zip(*) gets an empty list
				continue


			##############################################
			###		MULTI NUCLEOTIDE VARIANT CALLER    ###
			##############################################
			
			cposs = []  # Start positions of substitution groups
			crefs = []  # Merged reference bases
			cseqs = []  # Merged sequence bases
			clen = []   # Lengths of merged substitutions
			
			# Initialize the first group
			cposs.append(poss[0])
			ref_group = [refs[0]]
			seq_group = [seqs[0]]
			group_len = 1

			for i in range(1, len(poss)):
				if abs(poss[i] - poss[i - 1]) == 1: # check if neighbour to either side
					ref_group.append(refs[i])
					seq_group.append(seqs[i])
					group_len += 1
				else:
					# Store previous group
					crefs.append("".join(ref_group))
					cseqs.append("".join(seq_group))
					clen.append(group_len)
					# Start a new group
					cposs.append(poss[i])
					ref_group = [refs[i]]
					seq_group = [seqs[i]]
					group_len = 1
			
			# Add the last group
			crefs.append("".join(ref_group))
			cseqs.append("".join(seq_group))
			clen.append(group_len)
					

			for i in range(len(cposs)):
				end = cposs[i] + clen[i]
				mnv = f"{chrom}\t{cposs[i]}\t{end}\t{crefs[i]}\t{cseqs[i]}\t{clen[i]}"
				read_mnv = f"{name}\t{mnv}"

				# mnv_no_alt = f"{chrom}\t{cposs[i]}\t{end}\t{crefs[i]}"
				# no_alt_read_key = f"{name}\t{mnv_no_alt}" # does not hold alternative base
				# obs_base = cseqs[i]
				
				# all
				if not mnv in mnv_dict:
					mnv_dict[mnv] = 0
				mnv_dict[mnv] += 1

				# init dicts
				if not mnv in mnv_dict_consensus:
					mnv_dict_consensus[mnv] = 0
				# if not mnv_no_alt in mnv_dict_nonconsensus:
				# 	mnv_dict_nonconsensus[mnv_no_alt] = 1

				#non-consensus
				# if no_alt_read_key not in mnv_alt_dict: # track all obs_bases for position
				# 	mnv_alt_dict[no_alt_read_key] = set()
				# mnv_alt_dict[no_alt_read_key].add(obs_base)

				# consensus
				if not read_mnv in mnv_set:
					# first observation of mnv in read
					mnv_set.add(read_mnv)
				else:
					# second correct observation of mnv in read
					mnv_dict_consensus[mnv] += 1
					mnv_dict[mnv] -= 1 # consensus; hence adjust for counting twice
					mnv_set.discard(read_mnv)
					read_set.add(name)

				#check if non-consensus
				# if len(mnv_alt_dict[no_alt_read_key])>1:
				# 	mnv_dict_nonconsensus[mnv_no_alt] += 1
				# 	read_set_nonconsens.add(name)

			# if name in read_set:
			# 	if not fragment_length in frag_dict:
			# 		frag_dict[fragment_length] = 1
			# 	else:
			# 		frag_dict[fragment_length] += 1

	mt = 0
	ic = 0
	mc = 0
	nc = 0
	dc = 0

	# all
	for mnv, n_total in mnv_dict.items():
		n_consensus = mnv_dict_consensus.get(mnv, 0)
		n_nonconsensus = mnv_dict_nonconsensus.get(mnv, 0)
		# collect stats
		mt += n_total
		if n_consensus > 0:
			ic += 1
			mc += n_consensus
		
		if n_nonconsensus > 0:
			dc += 1
			nc += n_nonconsensus


		# format and write
		# print(f'{mnv}\t{n_consensus}\t{n_total}')
		record = f'{mnv}\t{n_consensus}\t{n_total}\n'
		fmnv.write(record.encode('utf-8'))

		# record_cons = f'{mnv}\t{n_consensus}\n'
		# fmnv.write(record_cons.encode('utf-8'))

		# record_noncons = f'{mnv}\t{n_nonconsensus}\n'
		# fmnv.write(record_noncons.encode('utf-8'))
		
	# consensus
	# for mnv, n in mnv_dict_consensus.items():
	# 	print(mnv, n)
	# 	if(n>0):
	# 		ic += 1
	# 		mc += n
	# 		record = f'{mnv}\t{n}\n'
	# 		# print(record)
	# 		fmnv.write(record.encode('utf-8'))

	# print summary stats
	#variants
	print("No. of read pairs with consensus MNVs:", len(read_set))
	print("No. of read pairs with non-consensus MNVs:", len(read_set_nonconsens))
	print("No. of MNVs (total):", len(mnv_dict))
	print("No. of MNVs (consensus):", ic)
	print("No. of MNVs (non-consensus):", dc)
	#mismatches
	print("No. of mismatches (total):", mt)
	print("No. of mismatches (consensus)", mc)
	print("Consensus/total mismatch ratio:", mc/mt)
	
	end_time = time.time()
	elapsed_time = end_time - start_time
	print("Saved consensus called MNVs to", output_mnv)
	print(f"Runtime: {elapsed_time:.4f} seconds")
	print("DONE")

	# for f, n in frag_dict.items():
	# 	print(f, n)
# print(mnv_alt_dict)
	# count 
	# awk 'FNR>1 {C[$7]++}END{for(i in C) print i, C[i]}' mnv_test.tsv

